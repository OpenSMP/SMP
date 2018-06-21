#include "SMP/SMPServer.hpp"
#include "SMP/network/net_io.hpp"
#include "SMP/Timer.hpp"
#include "SMP/literal.hpp"
#include "SMP/HElib.hpp"

#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/NumbTh.h>

//using boost::asio::ip::tcp;
static void randomize(Matrix &mat, long p = 3) {
    p = 11;
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = NTL::RandomBnd(p);
}

static long round_div(long a, long b) {
    return (a + b - 1) / b;
}

static void MyBuildLinPoly(std::vector<std::vector<zzX>> &Cs, const std::vector<NTL::ZZX> &Ls, const EncryptedArray *ea);

SMPServer::SMPServer() { }

SMPServer::~SMPServer() 
{
	if (context)
		delete context;
	if (ek)
		delete ek;
}

void SMPServer::run(tcp::iostream &conn, 
					const long n1,
					const long n2,
					const long n3) 
{
    setup(conn);
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);

	receive_ctx(conn);

    double evaluate_time = 0.;
    {
        AutoTimer timer(&evaluate_time);
        process_columns();
        evaluate();
    }
    extract_and_merge();
	response_ctx(conn);
    conn << evaluate_time;
}

void SMPServer::setup(tcp::iostream &conn)
{
	receive_context(conn, &context);
    ek = new FHEPubKey(*context);
    conn >> *ek;

    const EncryptedArray *ea = context->ea;
    const long d = ea->getDegree();
    const long l = ea->size();
    std::vector<NTL::ZZX> Ls(d);
    for (long i = 0; i < d; ++i)
        Ls[i] = NTL::ZZX(i, 1); // X^{d-1} -> X^i
    MyBuildLinPoly(linear_map_coeffs, Ls, context->ea);

    for (long i = 0; i < d; ++i) {
        size_t sze = linear_map_coeffs.at(i).size();
        std::vector<zzX> encoded_coeffs(sze);
        for (long j = 0; j < sze; ++j) {
            std::vector<zzX> tmp(l, linear_map_coeffs[i][j]);
            ea->encode(encoded_coeffs[j], tmp);
        }
        std::swap(linear_map_coeffs[i], encoded_coeffs);
    }
}

void SMPServer::process_columns()
{
    NTL::SetSeed(NTL::to_ZZ(123)); // use same seed for debugging
	randomize(A, ek->getPtxtSpace());
    randomize(B, ek->getPtxtSpace());
    ground_truth = mul(A, B);
    std::cout << ground_truth << "\n";

    const auto &ea = context->ea->getDerived(PA_zz_p());
    const long l = ea.size();
    const long d = ea.getDegree();
	const bool is_vec = A.NumRows() == 1;

	Matrix Bt;
	/// We compute A*B, but we use B tranpose.
	/// This allow us to write one internal::partition()
	/// for row-major case.
    transpose(&Bt, B);
    const long MAX_X2 = round_div(Bt.NumRows(), l);
    const long MAX_Y2 = round_div(Bt.NumCols(), d);
    const long MAX_Y1 = round_div(A.NumCols(), d);
	assert(MAX_Y1 == MAX_Y2);

    plain_B_blk.SetDims(MAX_X2, MAX_Y2);
    encoded_plain_B_blk.SetDims(MAX_X2, MAX_Y2);
    for (int y = 0; y < MAX_X2; y++) {
        for (int k = 0; k < MAX_Y2; k++) {
            internal::BlockId blk = {y, k};
            plain_B_blk[y][k] = internal::partition(Bt, blk, context->ea, true);
			if (is_vec) {
                ea.encode(encoded_plain_B_blk[y][k], plain_B_blk[y][k].polys);
			}
        }
    }
}

void SMPServer::receive_ctx(tcp::iostream &conn)
{
	const EncryptedArray *ea = context->ea;
    const long l = ea->size();
    const long d = ea->getDegree();
	const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);

    Ctxt tmp(*ek);
	enc_A_blk.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, tmp));
    for (int x = 0; x < MAX_X1; x++) {
        for (int k = 0; k < MAX_Y1; k++) {
            conn >> enc_A_blk[x][k];
        }
    }
}

void SMPServer::evaluate_mat_vec()
{
	const EncryptedArray *ea = context->ea;
	const long l = ea->size();
    const long d = ea->getDegree();
	assert(A.NumRows() == 1);
    const long MAX_Y1 = round_div(A.NumCols(), d);

	for (long col_B = 0; col_B < B.NumCols(); col_B += l) {
		long B_blk_idx = col_B / l;
		assert(B_blk_idx <= plain_B_blk.NumRows());
		Ctxt summation(*ek);
		for (long prtn = 0; prtn < MAX_Y1; prtn++) {
			Ctxt enc_blk(enc_A_blk.at(0).at(prtn));
			enc_blk.multByConstant(encoded_plain_B_blk[B_blk_idx][prtn]);
			summation += enc_blk;
		}
		results.push_back(summation);
	}
}

void SMPServer::evaluate()
{
    const auto& ea = context->ea->getDerived(PA_zz_p());
    const size_t l = static_cast<size_t>(ea.size());
    const size_t d = static_cast<size_t>(ea.getDegree());
	const size_t MAX_X1 = static_cast<size_t>(round_div(A.NumRows(), l));
    const size_t MAX_Y1 = static_cast<size_t>(round_div(A.NumCols(), d));
	const bool is_vec = A.NumRows() == 1;
	if (is_vec) {
		evaluate_mat_vec();
	} else {
        const size_t cols = B.NumCols();
		for (size_t A_blk_idx = 0; A_blk_idx < MAX_X1; A_blk_idx++) {
            const auto& enc_blk = enc_A_blk.at(A_blk_idx);
			for (long col_B = 0; col_B < cols; col_B++) {

				size_t B_blk_idx = col_B / l;
				size_t offset = col_B % l;
				assert(B_blk_idx <= plain_B_blk.NumRows());

                zzX plain_blk;
                Ctxt summation(*ek);
                const auto& row_blks = plain_B_blk[B_blk_idx];
				for (size_t prtn = 0; prtn < MAX_Y1; prtn++) {
					Ctxt eblk(enc_blk.at(prtn));
                    std::vector<zzX> _slots(l, row_blks[prtn].polys.at(offset));
                    ea.encode(plain_blk, _slots);
					eblk.multByConstant(plain_blk);
					summation += eblk;
				}
				results.push_back(summation);
			}
		}
	}
}

void SMPServer::extract_and_merge()
{
    std::list<Ctxt> merged_ctxts;
    const long batch_sze = context->ea->getDegree();
    long counter = 0;
    Ctxt ans(*ek);
    for (auto &ctx : results) {
        if (counter == batch_sze) {
            ans.modDownToLevel(1);
            merged_ctxts.push_back(ans);
            ans = Ctxt(*ek);
            counter = 0;
        }
        applyLinPolyLL(ctx, linear_map_coeffs.at(counter), batch_sze);
        ans += ctx;
        ++counter;
    }

    if (counter > 0) {
        ans.modDownToLevel(1);
        merged_ctxts.push_back(ans);
    }
    std::swap(results, merged_ctxts);
}

void SMPServer::response_ctx(tcp::iostream &conn)
{
	int64_t ctx_cnt = results.size();
	conn << ctx_cnt << std::endl;
	for (auto const& ctx : results)
		conn << ctx;
	conn.flush();
}

void MyBuildLinPoly(std::vector<std::vector<zzX>> &Cs,
                    const std::vector<NTL::ZZX> &Ls,
                    const EncryptedArray *ea)
{
    using RBak = PA_zz_p::RBak;
    using REbak = PA_zz_p::REBak;
    using RE = PA_zz_p::RE;
    using RX = PA_zz_p::RX;

    RBak bak; bak.save(); ea->restoreContext();
    REbak ebak; ebak.save(); ea->restoreContextForG();
    NTL::Lazy< NTL::Mat<RE> > linPolyMatrix;
    do {
        typename Lazy< Mat<RE> >::Builder builder(linPolyMatrix);
        if (!builder()) break;

        FHE_NTIMER_START(buildLinPolyCoeffs_invert);


        long p = ea->getPAlgebra().getP();
        long r = ea->getAlMod().getR();

        NTL::Mat<RE> M1;
        // build d x d matrix, d is taken from the current NTL context for G
        buildLinPolyMatrix(M1, p);
        NTL::Mat<RE> M2;
        ppInvert(M2, M1, p, r); // invert modulo prime-power p^r

        NTL::UniquePtr< Mat<RE> > ptr;
        ptr.make(M2);
        builder.move(ptr);
    } while (0);

    long d = ea->getDegree();
    Cs.resize(d);
    for (long i = 0; i < d; ++i) {
        RE L = NTL::conv<RE>(NTL::conv<RX>(Ls.at(i)));
        NTL::Vec<RE> tmp;
        NTL::mul(tmp, L, linPolyMatrix->operator[](d - 1));
        NTL::Vec<RX> _tmp;
        convert(_tmp, tmp);
        convert(Cs.at(i), _tmp);
    }
}
