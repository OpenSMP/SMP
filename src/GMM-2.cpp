/// General Matrix Multiplication
#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/NumbTh.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Matrix.hpp"
#include "CryptGMM/Timer.hpp"
#include "CryptGMM/HElib.hpp"
#include "CryptGMM/literal.hpp"
#include "CryptGMM/network/net_io.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <iostream>
#include <numeric>
#include <list>
using boost::asio::ip::tcp;
constexpr int REPEAT = 10;

inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

void zero(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = 0;
}

void randomize(Matrix &mat, long p = 3) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = NTL::RandomBnd(4);
}


void fill_compute(Matrix& mat,
				  long row_blk,
				  long col,
                  const std::vector<long> &inner_prod,
                  const EncryptedArray *ea)
{
    const long l = ea->size();
    assert(inner_prod.size() == l);
	const long row_start = row_blk * l;
    for (long ll = 0; ll < l; ll++) {
        long computed = inner_prod[ll];
        long row = row_start + ll;
		if (row < mat.NumRows()) {
			mat.put(row, col, computed);
		} else {
			break;
		}
    }
}

struct ClientBenchmark {
    std::vector<double> pack_times;
    std::vector<double> enc_times;
    std::vector<double> dec_times;
    std::vector<double> unpack_times;
    std::vector<double> total_times;
    int ctx_sent, ctx_recv;
};
ClientBenchmark clt_ben;

struct ServerBenchmark {
    std::vector<double> eval_times;
};
ServerBenchmark srv_ben;

void play_client(tcp::iostream &conn,
                 FHESecKey &sk,
                 FHEcontext &context,
                 const long n1,
                 const long n2,
                 const long n3) {
    FHEPubKey ek(sk);
	//* Convert to evalution key.
	//* This function is not provided by the origin HElib. Checkout our fork.
    ek.makeSymmetric();
    conn << ek;
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();

    NTL::SetSeed(NTL::to_ZZ(123));
    Matrix A, B, ground_truth;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A, ek.getPtxtSpace());
    randomize(B, ek.getPtxtSpace());
    ground_truth = mul(A, B);
    /// print grouth truth for debugging
    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
    const long MAX_X2 = round_div(B.NumCols(), l);

    std::vector<std::vector<Ctxt>> uploading;
    uploading.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, sk));
	double enc_time = 0.;
    double pack_time = 0.;
	/// encrypt matrix
	NTL::ZZX packed_poly;
	for (int x = 0; x < MAX_X1; x++) {
		for (int k = 0; k < MAX_Y1; k++) {
			internal::BlockId blk = {x, k};
			double one_pack_time, one_enc_time;
			auto block = internal::partition(A, blk, *ea, false);
			{/// packing
				AutoTimer timer(&one_pack_time);
				rawEncode(packed_poly, block.polys, context);
			}
			{/// encryption
				AutoTimer timer(&one_enc_time);
				sk.Encrypt(uploading[x][k], packed_poly);
			}
			pack_time += one_pack_time;
			enc_time += one_enc_time;
		}
	}
    clt_ben.pack_times.push_back(pack_time);
    clt_ben.enc_times.push_back(enc_time);

    /// send ciphertexts of matrix
    for (auto const& row : uploading) {
        for (auto const& ctx : row)
            conn << ctx;
    }
    clt_ben.ctx_sent = MAX_X1 * MAX_Y1;

    std::vector<GMMPrecompTable> tbls = precompute_gmm_tables(context);
    /// waiting results
    long rows_of_A = A.NumRows();
    long rows_of_Bt = B.NumCols(); // Bt::Rows = B::Cols
	int64_t ctx_cnt = 0;
	conn >> ctx_cnt;
    clt_ben.ctx_recv = ctx_cnt;
    std::vector<Ctxt> ret_ctxs(ctx_cnt, ek);
	for (size_t k = 0; k < ctx_cnt; k++)
		conn >> ret_ctxs.at(k);
    double eval_time = 0.;
    conn >> eval_time;
    srv_ben.eval_times.push_back(eval_time);
    /// decrypt
    Matrix computed;
    computed.SetDims(A.NumRows(), B.NumCols());
    zero(computed);
    int x = 0;
    int y = 0;
    std::vector<long> slots;
    std::vector<NTL::zz_pX> _slots;
    NTL::ZZX decrypted;
	double decrypt_time = 0.;
    double unpack_time = 0.;
	long ctx_idx = 0;
	bool dec_pass = true;
	for (const auto &ctx : ret_ctxs) {
		double one_dec_time, one_unpack_time;
		do {
			AutoTimer timer(&one_dec_time);
			dec_pass &= ctx.isCorrect();
			sk.Decrypt(decrypted, ctx);
		} while(0);
		do {
			AutoTimer timer(&one_unpack_time);
            extract_inner_products(slots, decrypted, tbls, context);
		} while(0);
        decrypt_time += one_dec_time;
        unpack_time += one_unpack_time;

		long row_blk = ctx_idx / B.NumCols();
		long column = ctx_idx % B.NumCols();
		ctx_idx += 1;
        fill_compute(computed, row_blk, column, slots, ea);
    }
    clt_ben.dec_times.push_back(decrypt_time);
    clt_ben.unpack_times.push_back(unpack_time);
	if (!::is_same(ground_truth, computed, NTL::zz_p::modulus()))
		std::cerr << "The computation seems wrong " << std::endl;
	if (!dec_pass)
		std::cerr << "Decryption might fail" << std::endl;
}

void play_server(tcp::iostream &conn,
                 const long n1,
                 const long n2,
                 const long n3) {
    auto context = receive_context(conn);
    NTL::zz_p::init(context.zMStar.getP());
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    FHEPubKey ek(context);
    conn >> ek;

    NTL::SetSeed(NTL::to_ZZ(123)); /// use same seed for debugging
    Matrix A, B, ground_truth;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A, ek.getPtxtSpace());
    randomize(B, ek.getPtxtSpace());
    ground_truth = mul(A, B);
    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);

    Matrix Bt;
	/// We compute A*B, but we use B tranpose.
	/// This allow us to write one internal::partition()
	/// for row-major case.
    transpose(&Bt, B);
    const long MAX_X2 = round_div(Bt.NumRows(), l);
    const long MAX_Y2 = round_div(Bt.NumCols(), d);
	assert(MAX_Y1 == MAX_Y2);
    NTL::Mat<internal::PackedRows> plain_B_blk; // 2D-array of polynomials
    plain_B_blk.SetDims(MAX_X2, MAX_Y2);
    for (int y = 0; y < MAX_X2; y++) {
        for (int k = 0; k < MAX_Y2; k++) {
            internal::BlockId blk = {y, k};
            plain_B_blk[y][k] = internal::partition(Bt, blk, *ea, true);
        }
    }

    /// receving ciphertexts from the client
    std::vector<std::vector<Ctxt>> enc_A_blk;
    enc_A_blk.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, ek));
    for (int x = 0; x < MAX_X1; x++) {
        for (int k = 0; k < MAX_Y1; k++)
            conn >> enc_A_blk[x][k];
    }

    /// compute the matrix mulitplication
    double computation{0.};
	std::list<Ctxt> results;
	do {
		AutoTimer timer(&computation);
		for (long A_blk_idx = 0; A_blk_idx < MAX_X1; A_blk_idx++) {
			for (long col_B = 0; col_B < B.NumCols(); col_B++) {
				long B_blk_idx = col_B / l;
				long offset = col_B % l;
				assert(B_blk_idx <= plain_B_blk.NumRows());
				Ctxt summation(ek);
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
	} while (0);

	int64_t ctx_cnt = results.size();
	conn << ctx_cnt << std::endl;
	for (auto const& ctx : results)
		conn << ctx;
    /// sent the evalution time, just for statistics
    conn << computation;
}

int run_client(std::string const& addr, long port,
               long n1, long n2, long n3) {
    const long m = 8192;
    const long p = 84737;
    const long r = 1;
    const long L = 2;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 51;
    buildModChain(context, L);
    std::cerr << "kappa = " << context.securityLevel() << std::endl;
    std::cerr << "slot = " << context.ea->size() << std::endl;
    std::cerr << "degree = " << context.ea->getDegree() << std::endl;
    FHESecKey sk(context);
    sk.GenSecKey(64);

    for (long t = 0; t < REPEAT; t++) {
        tcp::iostream conn(addr, std::to_string(port));
        if (!conn) {
            std::cerr << "Can not connect to server!" << std::endl;
            return -1;
        }

        /// send FHEcontext obj
        double all_time;
        do {
            send_context(conn, context);
            AutoTimer time(&all_time);
            /// send the evaluation key
            play_client(conn, sk, context, n1, n2, n3);
        } while(0);
        clt_ben.total_times.push_back(all_time);
        conn.close();
    }
    return 1;
}

int run_server(long port, long n1, long n2, long n3) {
    boost::asio::io_service ios;
    tcp::endpoint endpoint(tcp::v4(), port);
    tcp::acceptor acceptor(ios, endpoint);

    for (long run = 0; run < REPEAT; run++) {
        tcp::iostream conn;
        boost::system::error_code err;
        acceptor.accept(*conn.rdbuf(), err);
        if (!err) {
            play_server(conn, n1, n2, n3);
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    ArgMapping argmap;
    long role = -1;
    long n1 = 8;
    long n2 = 8;
    long n3 = 8;
    std::string addr = "127.0.0.1";
	long port = 12345;
    argmap.arg("N", n1, "n1");
    argmap.arg("M", n2, "n2");
    argmap.arg("D", n3, "n3");
    argmap.arg("R", role, "role. 0 for server and 1 for client");
	argmap.arg("a", addr, "server address");
	argmap.arg("p", port, "port");
    argmap.parse(argc, argv);
    if (role == 0) {
        run_server(port, n1, n2, n3);
    } else if (role == 1) {
        run_client(addr, port, n1, n2, n3);
        auto time = mean_std(clt_ben.pack_times);
        printf("%.3f %.3f ", time.first, time.second);

        time = mean_std(clt_ben.enc_times);
        printf("%.3f %.3f ", time.first, time.second);

        time = mean_std(clt_ben.dec_times);
        printf("%.3f %.3f ", time.first, time.second);

        time = mean_std(clt_ben.unpack_times);
        printf("%.3f %.3f ", time.first, time.second);

        time = mean_std(clt_ben.total_times);
        printf("%.3f %.3f ", time.first, time.second);

        time = mean_std(srv_ben.eval_times);
        printf("%.3f %.3f ", time.first, time.second);
        printf("%d %d\n", clt_ben.ctx_sent, clt_ben.ctx_recv);
    } else {
		argmap.usage("General Matrix Multiplication for |N*M| * |M*D|");
		return -1;
	}
}
