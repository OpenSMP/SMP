#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/NumbTh.h>
#include <NTL/lzz_pX.h>

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
static void randomize(Matrix &mat);
static void randomize(NTL::ZZX &poly, long deg, long p);

void pack_rows(std::vector<NTL::ZZX> &polys, Matrix const& mat, const long m)
{
	const long phim = phi_N(m);
	assert(phim == (m >> 1));
	assert(phim >= mat.NumCols());
	for (long row = 0; row < mat.NumRows(); row++) {
		NTL::ZZX poly;
		poly.SetLength(phim);
		auto dst = poly.rep.begin();
		auto itr = mat[row].begin();
		auto end = mat[row].end();
		while (itr != end)
			*dst++ = NTL::to_ZZ(*itr++);
		polys.push_back(poly);
	}
}

void pack_column(NTL::ZZX &poly, Matrix const& mat, const long col, const long m)
{
	const long phim = phi_N(m);
	assert(phim == (m >> 1));
	assert(phim >= mat.NumRows());
	assert(col >= 0 && col < mat.NumCols());
	poly.SetLength(phim);
	for (long r = 0; r < mat.NumRows(); r++) {
		long idx = (phim - r) % phim;
		long sig = r > 0 ? -1 : 1;
		if (r > 0)
			NTL::SetCoeff(poly, idx, -mat[r][col]);
		else
			NTL::SetCoeff(poly, idx, mat[r][col]);
	}
}

std::vector<long> get_inner_products(NTL::ZZX const& poly, const long vec_len) {
	std::vector<long> ips;
	for (long i = 0; i <= NTL::deg(poly); i += vec_len)
		ips.push_back(NTL::to_long(poly[i]));
	return ips;
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

	std::vector<NTL::ZZX> polys;
	double pack_time;
	do {
		AutoTimer timer(&pack_time);
		pack_rows(polys, A, context.zMStar.getM());
	} while (0);
	clt_ben.pack_times.push_back(pack_time);

	std::vector<Ctxt> ctxts;
	double enc_time;
	do {
		AutoTimer timer(&enc_time);
		std::transform(polys.cbegin(), polys.cend(), std::back_inserter(ctxts),
					   [&sk](NTL::ZZX const& poly) -> Ctxt {
					   Ctxt ctx(sk);
					   sk.Encrypt(ctx, poly);
					   return ctx;
					   });
	} while (0);
	clt_ben.enc_times.push_back(enc_time);

    conn << static_cast<int64_t>(ctxts.size()) << std::endl;
	for (const auto& ctx : ctxts)
		conn << ctx;
    clt_ben.ctx_sent = ctxts.size();

	int64_t result_ctx_cnt;
	conn >> result_ctx_cnt;
    clt_ben.ctx_recv = result_ctx_cnt;
    std::vector<Ctxt> results(result_ctx_cnt, ek);
	for (size_t i = 0; i < result_ctx_cnt; i++)
		conn >> results[i];
    double eval_time;
    conn >> eval_time;
    srv_ben.eval_times.push_back(eval_time);

	double dec_time;
	std::vector<NTL::ZZX> decrypted;
	do {
		AutoTimer timer(&dec_time);
		std::transform(results.cbegin(), results.cend(), std::back_inserter(decrypted),
					   [&sk](Ctxt const& ctx) -> NTL::ZZX {
					   NTL::ZZX dec;
					   if (!ctx.isCorrect())
					       std::cout << "decryption might fail" << std::endl;
					   sk.Decrypt(dec, ctx);
					   return dec;
					   });
	} while (0);
	clt_ben.dec_times.push_back(dec_time);

	double unpack_time;
	do {
		AutoTimer timer(&unpack_time);
		for (auto &dec : decrypted)
			auto ips = get_inner_products(dec, B.NumRows());
	} while(0);
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
	std::vector<Ctxt> ctxts(ctx_cnt, evk);
	for (size_t i = 0; i < ctx_cnt; i++) // receive ciphertexts from the client.
		conn >> ctxts[i];

	double eval_time;
	std::list<Ctxt> results;
	do {
		AutoTimer timer(&eval_time);
		for (long col = 0; col < B.NumCols(); col++) {
			NTL::ZZX poly;
			pack_column(poly, B, col, context.zMStar.getM());
			auto tmp(ctxts);
			for (auto &ctx : tmp) {
				ctx.multByConstant(poly);
				NTL::ZZX rnd; // adding the random share
				randomize(rnd, context.zMStar.getM(), context.zMStar.getP());
				ctx.addConstant(rnd);
				results.push_back(ctx);
			}
		}
	} while (0);
	srv_ben.eval_times.push_back(eval_time);

    conn << static_cast<int64_t>(results.size()) << std::endl;
    for (const auto &ctx : results) // sending result back
		conn << ctx;
}

int run_client(std::string const& addr, long port,
			   long n1, long n2, long n3) {
    const long m = 8192;
    const long p = 769;
    const long r = 2;
    const long L = 2;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 30 + std::ceil(std::log(m)/2 + r * std::log(p));

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
		run_client(addr, port, n1, n2, n3);

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
	save_matrix(std::cout, AB);

	std::vector<NTL::ZZX> polys;
	pack_rows(polys, A, context.zMStar.getM());
	std::vector<Ctxt> ctxts;
	std::transform(polys.cbegin(), polys.cend(), std::back_inserter(ctxts),
				   [&sk](NTL::ZZX const& poly) -> Ctxt {
					   Ctxt ctx(sk);
					   sk.Encrypt(ctx, poly);
					   return ctx;
				   });

	std::list<Ctxt> results;
	for (long col = 0; col < B.NumCols(); col++) {
		NTL::ZZX poly;
		pack_column(poly, B, col, context.zMStar.getM());
		auto tmp(ctxts);
		for (auto &ctx : tmp) {
			ctx.multByConstant(poly);
			NTL::ZZX rnd; // adding the random share
			randomize(rnd, context.zMStar.getM(), context.zMStar.getP());
			ctx.addConstant(rnd);
			results.push_back(ctx);
		}
	}

	std::vector<NTL::ZZX> decrypted;
	std::transform(results.cbegin(), results.cend(), std::back_inserter(decrypted),
				   [&sk](Ctxt const& ctx) -> NTL::ZZX {
				   NTL::ZZX dec;
				   sk.Decrypt(dec, ctx);
				   return dec;
				   });
	for (auto &dec : decrypted) {
		auto ips = get_inner_products(dec, B.NumRows());
	}
	return 0;
}
#endif
void randomize(Matrix &mat) {
	for (long i = 0; i < mat.NumRows(); i++)
		for (long j = 0; j < mat.NumCols(); j++)
			mat[i][j] = NTL::RandomBnd(10L);
}

void randomize(NTL::ZZX &poly, long m, long p) {
	NTL::zz_p::init(p);
	NTL::zz_pX ply;
	long phim = phi_N(m);
	assert(phim == (m >> 1));
	NTL::random(ply, phim);
	NTL::conv(poly, ply);
}

