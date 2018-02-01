#include "dgk.h"
#include "SMP/Matrix.hpp"
#include "SMP/Timer.hpp"
#include "HElib/NumbTh.h"
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <iostream>
#include <numeric>
#include <list>
using boost::asio::ip::tcp;
constexpr int REPEAT = 1;

void get_rand_file( void* buf, int len, char* file )
{
	FILE* fp;
	void* p;

	fp = fopen(file, "r");

	p = buf;
	while( len )
	{
		size_t s;
		s = fread(p, 1, len, fp);
		p += s;
		len -= s;
	}

	fclose(fp);
}

void get_rand_devrandom( void* buf, int len )
{
	get_rand_file(buf, len, "/dev/random");
}

void get_rand_devurandom( void* buf, int len )
{
	get_rand_file(buf, len, "/dev/urandom");
}

using get_rand_t = void (*)( void* buf, int len );

void init_rand( gmp_randstate_t rand, get_rand_t get_rand, int bytes )
{
	void* buf;
	mpz_t s;

	buf = malloc(bytes);
	get_rand(buf, bytes);

	gmp_randinit_default(rand);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	gmp_randseed(rand, s);
	mpz_clear(s);

	free(buf);
}

using EncMatrix = NTL::Mat<mpz_t>;

void randomize_matrix(Matrix &mat)
{
	constexpr long domain = (1 << 16);
	for (long r = 0; r < mat.NumRows(); r++) {
		for (long c = 0; c < mat.NumCols(); c++)
			mat[r][c] = NTL::RandomBnd(domain);
	}
}

void encrypt_matrix(EncMatrix &dst, Matrix const& src, dgk_pubkey_t const* pk)
{
	if (src.NumRows() != dst.NumRows() or src.NumCols() != dst.NumCols())
		return;
	gmp_randstate_t gmp_rand;
    init_rand(gmp_rand, get_rand_devurandom, pk->bits / 8 + 1);
	mpz_t pt;
	mpz_init(pt);
	for (long r = 0; r < src.NumRows(); r++) {
		for (long c = 0; c < src.NumCols(); c++) {
			mpz_set_ui(pt, src[r][c]);
			mpz_init(dst[r][c]);
			dgk_encrypt_plain(dst[r][c], pk, pt, gmp_rand);
		}
	}
	mpz_clear(pt);
}

struct ClientBenchmark {
    std::vector<double> enc_times;
    std::vector<double> dec_times;
    std::vector<double> total_times;
    int ctx_sent, ctx_recv;
};
ClientBenchmark clt_ben;

struct ServerBenchmark {
    std::vector<double> eval_times;
};
ServerBenchmark srv_ben;


void play_server(std::iostream &conn, long n1, long n2, long n3);
void play_client(std::iostream &conn, dgk_prvkey_t *sk, dgk_pubkey_t *pk, long n1, long n2, long n3);

void run_server(long port, long n1, long n2, long n3)
{
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
}

int run_client(std::string const& addr, long port, long n1, long n2, long n3) {
	dgk_prvkey_t *sk;
	dgk_pubkey_t *pk;
	dgk_keygen(1024, 16, &pk, &sk); // 16-bit plaintexts
    for (long t = 0; t < REPEAT; t++) {
        tcp::iostream conn(addr, std::to_string(port));
        if (!conn) {
            std::cerr << "Can not connect to server!" << std::endl;
            return -1;
        }

        double all_time = 0.;
        do {
            AutoTimer time(&all_time);
            play_client(conn, sk, pk, n1, n2, n3);
        } while(0);
        clt_ben.total_times.push_back(all_time);
        conn.close();
    }
	dgk_freeprvkey(sk);
	dgk_freepubkey(pk);
    return 1;
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
    } else {
		argmap.usage("Secure Matrix Multiplication for |N*M| * |M*D|");
		return -1;
	}
	
	return 0;
}

void play_server(std::iostream &conn, long n1, long n2, long n3)
{
	Matrix A, B;
	NTL::SetSeed(NTL::to_ZZ(123)); // for debugging
	A.SetDims(n1, n2);
	B.SetDims(n2, n3);
	randomize_matrix(A);
	randomize_matrix(B);


	std::vector<uint32_t> buffer(2048 >> 5);
	dgk_pubkey_t pk;
	receive_pk(&pk, buffer, conn);
	
	gmp_randstate_t gmp_rand;
    init_rand(gmp_rand, get_rand_devurandom, pk.bits / 8 + 1);
    
 	mpz_t ct, pt;
	mpz_init(ct);
	mpz_init(pt);
	mpz_set_ui(pt, 10);
	dgk_encrypt_plain(ct, &pk, pt, gmp_rand);
	send_mpz(ct, buffer, conn);
}

void play_client(std::iostream &conn, dgk_prvkey_t *sk, dgk_pubkey_t *pk, long n1, long n2, long n3)
{
	std::vector<uint32_t> buffer(2048 >> 5);
	send_pk(pk, buffer, conn);
	mpz_t ct;
	mpz_init(ct);
	receive_mpz(ct, buffer, conn);
	mpz_t pt;
	mpz_init(pt);
	dgk_decrypt(pt, pk, sk, ct);
	gmp_printf("%Zd\n", pt);
}
