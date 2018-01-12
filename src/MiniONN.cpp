#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/NumbTh.h>

#include "CryptGMM/Matrix.hpp"
#include "CryptGMM/Timer.hpp"
#include "CryptGMM/literal.hpp"
#include "CryptGMM/network/net_io.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <iostream>
#include <numeric>
#include <list>
using boost::asio::ip::tcp;
constexpr long REPEATS = 10;

struct EncVec {
	long length;
	std::vector<Ctxt> ctxts;
};

void pack_vector(std::vector<NTL::ZZX> &out, Matrix const& mat, long COL, EncryptedArray const* ea)
{
	const long ROWS = mat.NumRows();
	long num_polys = (ROWS + ea->size() - 1) / ea->size();
	long itr = 0;
	out.resize(num_polys);
	for (long i = 0; i < num_polys; i++) {
		std::vector<long> slots(ea->size(), 0);
		for (long j = 0; j < ea->size() and itr < ROWS; j++) {
			slots[j] = mat[itr++][COL];
		}
		ea->encode(out[i], slots);
	}
}

void rotate_vector(Matrix &vec, long offset)
{
	assert(vec.NumCols() == 1);
	long length = vec.NumRows();
	while (offset < 0)
		offset += length;
	if (offset == 0)
		return;
	Matrix tmp;
	tmp.SetDims(vec.NumRows(), 1);
	for (long i = 0; i < length; i++)
		tmp[i][0] = vec[(i + offset) % length][0];
	vec = tmp;
}

void encrypt_vector(EncVec &out,
					Matrix const& vec,
					FHESecKey const& key,
					double *pack_time = nullptr)
{
	assert(vec.NumCols() == 1);
	FHEcontext const& context= key.getContext();
	EncryptedArray const* ea = context.ea;
	out.length = vec.NumRows();
	std::vector<NTL::ZZX> polys;
	{
		AutoTimer timer(pack_time);
		pack_vector(polys, vec, 0, ea);
	}
	out.ctxts.resize(polys.size(), key);
	for (size_t i = 0; i < polys.size(); i++) {
		key.Encrypt(out.ctxts[i], polys.at(i));
	}
}

void columns_to_vector(Matrix &vec, Matrix const&mat, long from, long to)
{
	vec.SetDims((to - from) * mat.NumRows(), 1);
	for (long j = from; j < to; j++) {
		long offset = (j - from) * mat.NumRows();
		for (long i = 0; i < mat.NumRows(); i++) {
			vec[offset + i][0] = mat[i][j];
		}
	}
}

void rows_to_vector(Matrix &vec, Matrix const&mat, long from, long to)
{
	assert(from >= 0 && to > from && to <= mat.NumRows());
	const long length = (to - from) * mat.NumCols();
	vec.SetDims(length, 1);
	for (long i = from; i < to; i++) {
		long offset = (i - from) * mat.NumCols();
		for (long j = 0; j < mat.NumCols(); j++) {
			vec[offset + j][0] = mat[i][j];
		}
	}
}

/// To fully pack the matrix, we concat each row of matrix into a long vector.
void encrypt_matrix_as_vector(EncVec &out,
							  Matrix const& mat,
							  FHESecKey const& key,
							  double *pack_time = nullptr)
{
	Matrix vec;
	rows_to_vector(vec, mat, 0, mat.NumRows());
	encrypt_vector(out, vec, key, pack_time);
}

void randomize(Matrix &mat) {
	for (long i = 0; i < mat.NumRows(); i++)
		for (long j = 0; j < mat.NumCols(); j++)
			mat[i][j] = NTL::RandomBnd(10L);
}

void mat_mult(std::list<Ctxt> &out,
			  EncVec const& enc_mat,
			  Matrix const& mat,
			  EncryptedArray const *ea)
{
	const auto row_cnt = mat.NumRows();
	const auto col_cnt = mat.NumCols();
	long cols_per_ctxt = ea->size() / row_cnt;
	assert(cols_per_ctxt >= 1 && "Current version only work for n2 <= ea->size()");
	for (long col = 0; col < col_cnt; col += cols_per_ctxt) {
		long cols_to_pack = std::min(col_cnt - col, cols_per_ctxt);
		Matrix vec;
		columns_to_vector(vec, mat, col, col + cols_to_pack);

		Matrix rnd;
		rnd.SetDims(cols_to_pack * row_cnt, 1);
		for (long rot = 0; rot < cols_to_pack; rot++) { // rotatiton
			long offset = rot * row_cnt;
			rotate_vector(vec, offset);
			std::vector<NTL::ZZX> polys;
			pack_vector(polys, vec, 0, ea); // pack the 0-th column
			assert(polys.size() == 1);

			for (const auto &ctx : enc_mat.ctxts) {
				Ctxt tmp(ctx);
				tmp.multByConstant(polys[0]); // multiplication

                randomize(rnd);
				std::vector<NTL::ZZX> rnd_poly;
				pack_vector(rnd_poly, rnd, 0, ea);
				assert(rnd_poly.size() == 1);
                tmp.addConstant(rnd_poly[0]); // adding the random share

				out.push_back(std::move(tmp));
			}
		}
	}
}

struct ClientBenchmark {
    std::vector<double> pack_times;
	std::vector<double> enc_times;
    std::vector<double> unpack_times;
	std::vector<double> dec_times;
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
                 const long n3)
{
	FHEPubKey ek(sk);
    ek.makeSymmetric();
    conn << ek; // send evaluation key
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();

    NTL::SetSeed(NTL::to_ZZ(123)); //use fixed seed for debugging
    Matrix A, B;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A);
    randomize(B);

	EncVec enc_vec;
	double pack_time;
	double enc_time;
	do {
		AutoTimer timer(&enc_time);
		encrypt_matrix_as_vector(enc_vec, A, sk, &pack_time);
	} while(0);
	clt_ben.pack_times.push_back(pack_time);
	clt_ben.enc_times.push_back(enc_time);

    conn << static_cast<int64_t>(enc_vec.ctxts.size()) << std::endl;
    for (const auto &ctx : enc_vec.ctxts)
		conn << ctx;
    clt_ben.ctx_sent = enc_vec.ctxts.size();

	int64_t result_ctx_cnt;
	conn >> result_ctx_cnt;
    clt_ben.ctx_recv = result_ctx_cnt;
    std::vector<Ctxt> result(result_ctx_cnt, ek);
	for (size_t i = 0; i < result_ctx_cnt; i++)
		conn >> result[i];
    double eval_time = 0.;
    conn >> eval_time;
    srv_ben.eval_times.push_back(eval_time);

	double dec_time;
	double unpack_time;
	std::vector<long> slots;
	NTL::ZZX decrypted;
	for (const auto &ctx : result) {
		double one_dec_time;
		double one_unpack_time;
		{
			AutoTimer timer(&one_dec_time);
            if (!ctx.isCorrect())
				std::cerr << "decryption might fail" << std::endl;
			sk.Decrypt(decrypted, ctx);
		}
		{
			AutoTimer timer(&one_unpack_time);
            ea->decode(slots, decrypted);
		}
		dec_time += one_dec_time;
		unpack_time += one_unpack_time;
	}
	clt_ben.dec_times.push_back(dec_time);
	clt_ben.unpack_times.push_back(unpack_time);
}

void play_server(tcp::iostream &conn,
                 const long n1,
                 const long n2,
                 const long n3) {
    FHEcontext context = receive_context(conn);
	FHEPubKey evk(context);
	conn >> evk;

	NTL::SetSeed(NTL::to_ZZ(123)); //use fixed seed for debugging

    Matrix A, B;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A);
    randomize(B);

	int64_t ctx_cnt;
	conn >> ctx_cnt;
	EncVec enc_vec;
	enc_vec.ctxts.resize(ctx_cnt, evk);
	for (size_t i = 0; i < ctx_cnt; i++) // receive ciphertexts from the client.
		conn >> enc_vec.ctxts[i];

	double eval_time;
	std::list<Ctxt> result;
	{
		AutoTimer timer(&eval_time);
		mat_mult(result, enc_vec, B, context.ea);
	}
	srv_ben.eval_times.push_back(eval_time);

    conn << static_cast<int64_t>(result.size()) << std::endl;
    for (const auto &ctx : result) // sending result back
		conn << ctx;
    conn << eval_time; // send back eval time for statistics
}

int run_client(std::string const& addr, long port,
			   long n1, long n2, long n3) {
    const long m = 8192;
    const long p = 65537;
    const long r = 1;
    const long L = 3;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 20; // trial and error to find this value.
    buildModChain(context, L);
    std::cerr << "kappa = " << context.securityLevel() << std::endl;
    std::cerr << "slot = " << context.ea->size() << std::endl;
    std::cerr << "degree = " << context.ea->getDegree() << std::endl;
    FHESecKey sk(context);
    sk.GenSecKey(64);
    for (long t = 0; t < REPEATS; t++) {
        tcp::iostream conn(addr, std::to_string(port));
        if (!conn) {
            std::cerr << "Can not connect to server!" << std::endl;
            return -1;
        }
        /// send FHEcontext obj
        send_context(conn, context);
        double all_time;
        do {
            AutoTimer time(&all_time);
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

    for (long run = 0; run < REPEATS; run++) {
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
        int st = run_client(addr, port, n1, n2, n3);
		auto times = mean_std(clt_ben.pack_times);
		printf("%.3f %.3f ", times.first, times.second);

		times = mean_std(clt_ben.enc_times);
		printf("%.3f %.3f ", times.first, times.second);

        times = mean_std(clt_ben.dec_times);
        printf("%.3f %.3f ", times.first, times.second);

		times = mean_std(clt_ben.unpack_times);
		printf("%.3f %.3f ", times.first, times.second);

        times = mean_std(clt_ben.total_times);
        printf("%.3f %.3f ", times.first, times.second);

        times = mean_std(srv_ben.eval_times);
        printf("%.3f %.3f ", times.first, times.second);
        printf("%d %d\n", clt_ben.ctx_sent, clt_ben.ctx_recv);
    } else {
		argmap.usage("General Matrix Multiplication for |N*M| * |M*D|");
		return -1;
	}
}

#if 0
int main(int argc, char *argv[]) {
	ArgMapping argmap;
    long n1 = 8;
    long n2 = 8;
	long n3 = 8;
    argmap.arg("N", n1, "n1");
    argmap.arg("M", n2, "n2");
    argmap.arg("K", n3, "n3");
    argmap.parse(argc, argv);

	Matrix A, B;
	A.SetDims(n1, n2);
	randomize(A);
	B.SetDims(n2, n3);
	randomize(B);
	Matrix AB = mul(A, B);

	FHEcontext context(8192, 8191, 1);
	std::cout << "slots " << context.ea->size() << std::endl;
	buildModChain(context, 4);
	std::cout << "security level " << context.securityLevel() << std::endl;
	FHESecKey sk(context);
	sk.GenSecKey(64);
	auto ea = context.ea;
	EncVec enc_vec;
	double pack_time;
	double enc_time;
	{
		AutoTimer timer(&enc_time);
		encrypt_matrix_as_vector(enc_vec, A, sk, &pack_time);
	}

	std::list<Ctxt> result;
	double eval_time;
	{
		AutoTimer timer(&eval_time);
		mat_mult(result, enc_vec, B, ea);
	}
	double dec_time;
	{
		AutoTimer timer(&dec_time);
		std::vector<long> slots;
		for (auto &ctx : result) {
			ea->decrypt(ctx, sk, slots);
		}
	}
	printf("%.3f %.3f %.3f %zd %zd\n",
		   enc_time, dec_time, eval_time,
		   enc_vec.ctxts.size(), result.size());
	return 0;
}
#endif
