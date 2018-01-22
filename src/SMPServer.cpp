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
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = NTL::RandomBnd(p);
}

static long round_div(long a, long b) {
    return (a + b - 1) / b;
}

std::vector<double> SMPServer::setup_times;
std::vector<double> SMPServer::process_columns_times;
std::vector<double> SMPServer::receive_ctx_times;
std::vector<double> SMPServer::evaluate_times;
std::vector<double> SMPServer::response_ctx_times;

SMPServer::~SMPServer() 
{
	if (context)
		delete context;
	if (ek)
		delete ek;
}

void SMPServer::print_statistics() 
{
	double total = 0.;
	printf("setup process_columns receive_ctx evaluate response_ctx total\n");
	auto time = mean_std(setup_times);
	printf("%.3f ", time.first);
	total += time.first;

	time = mean_std(process_columns_times);
	printf("%.3f ", time.first);
	total += time.first;

	time = mean_std(receive_ctx_times);
	printf("%.3f ", time.first);
	total += time.first;

	time = mean_std(evaluate_times);
	printf("%.3f ", time.first);
	total += time.first;

	time = mean_std(response_ctx_times);
	printf("%.3f ", time.first);
	total += time.first;
	
	printf(": %.3f\n", total);
}

void SMPServer::run(tcp::iostream &conn, 
					const long n1,
					const long n2,
					const long n3) 
{
	A.SetDims(n1, n2);
	B.SetDims(n2, n3);
	setup(conn);
	process_columns();
	receive_ctx(conn);
	evaluate();
	response_ctx(conn);
}

void SMPServer::setup(tcp::iostream &conn)
{
	setup_times.push_back(0.);
	AutoTimer timer(&(setup_times.back()));
	receive_context(conn, &context);
	NTL::zz_p::init(context->zMStar.getP());
    ek = new FHEPubKey(*context);
    conn >> *ek;

	NTL::SetSeed(NTL::to_ZZ(123)); // use same seed for debugging
	randomize(A, ek->getPtxtSpace());
    randomize(B, ek->getPtxtSpace());
    ground_truth = mul(A, B);
}

void SMPServer::process_columns()
{
	process_columns_times.push_back(0.);
	AutoTimer timer(&(process_columns_times.back()));
	const EncryptedArray *ea = context->ea;
    const long l = ea->size();
    const long d = ea->getDegree();
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
            plain_B_blk[y][k] = internal::partition(Bt, blk, ea, true);
			if (is_vec) {
				// NTL::zz_pX encoded;
				rawEncode(encoded_plain_B_blk[y][k],
						  plain_B_blk[y][k].polys, *context);
				// NTL::conv(, encoded);
			}
        }
    }
}

void SMPServer::receive_ctx(tcp::iostream &conn)
{
	receive_ctx_times.push_back(0.);
	AutoTimer timer(&(receive_ctx_times.back()));
	const EncryptedArray *ea = context->ea;
    const long l = ea->size();
    const long d = ea->getDegree();
	const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);

	enc_A_blk.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, *ek));
    for (int x = 0; x < MAX_X1; x++) {
        for (int k = 0; k < MAX_Y1; k++)
            conn >> enc_A_blk[x][k];
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
		summation.modDownToLevel(1);
		results.push_back(summation);
	}
}

void SMPServer::evaluate()
{
	evaluate_times.push_back(getTimerByName("FROM_POLY_OUTPUT")->getTime() * 1000.);
	AutoTimer timer(&(evaluate_times.back())); // measure evaluation time
	const EncryptedArray *ea = context->ea;
    const long l = ea->size();
    const long d = ea->getDegree();
	const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
	const bool is_vec = A.NumRows() == 1;
	if (is_vec) {
		evaluate_mat_vec();
	} else {
		for (long A_blk_idx = 0; A_blk_idx < MAX_X1; A_blk_idx++) {
			for (long col_B = 0; col_B < B.NumCols(); col_B++) {
				long B_blk_idx = col_B / l;
				long offset = col_B % l;
				assert(B_blk_idx <= plain_B_blk.NumRows());
				Ctxt summation(*ek);
				for (long prtn = 0; prtn < MAX_Y1; prtn++) {
					Ctxt enc_blk(enc_A_blk.at(A_blk_idx).at(prtn));
					NTL::ZZX plain_blk;
					NTL::conv(plain_blk, plain_B_blk[B_blk_idx][prtn].polys.at(offset));
					enc_blk.multByConstant(plain_blk);
					summation += enc_blk;
				}
				summation.modDownToLevel(1);
				results.push_back(summation);
			}
		}
	}
}

void SMPServer::response_ctx(tcp::iostream &conn)
{
	response_ctx_times.push_back(0.);
	AutoTimer timer(&(response_ctx_times.back()));
	int64_t ctx_cnt = results.size();
	conn << ctx_cnt << std::endl;
	for (auto const& ctx : results)
		conn << ctx;
    /// sent the evalution time, just for statistics
	evaluate_times.back() += getTimerByName("TO_POLY_OUTPUT")->getTime() * 1000.;
    conn << evaluate_times.back() + process_columns_times.back();
	conn.flush();
}
