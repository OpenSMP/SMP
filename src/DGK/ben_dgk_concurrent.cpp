#include "dgk.h"
#include "SMP/Matrix.hpp"
#include "SMP/Timer.hpp"
#include "SMP/network/ConcurrentQueue.hpp"
#include "HElib/NumbTh.h"
#include "SMP/literal.hpp"
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <iostream>
#include <numeric>
#include <list>
#include <thread>
using boost::asio::ip::tcp;

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

using EncMatrix = std::vector<mpz_t>;

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
	int rows = src.NumRows();
	int cols = src.NumCols();
	if (dst.size() != rows * cols)
		return;
	gmp_randstate_t gmp_rand;
    init_rand(gmp_rand, get_rand_devurandom, pk->bits / 8 + 1);
	mpz_t pt;
	mpz_init(pt);
	for (long r = 0; r < rows; r++) {
		int offset = r * cols;
		for (long c = 0; c < cols; c++) {
			mpz_set_ui(pt, src[r][c]);
			mpz_init(dst[offset + c]);
			dgk_encrypt_plain(dst[offset + c], pk, pt, gmp_rand);
		}
	}
	mpz_clear(pt);
}

void play_server(std::iostream &conn, long n1, long n2, long n3);
void play_client(std::iostream &conn, dgk_prvkey_t *sk, dgk_pubkey_t *pk, long n1, long n2, long n3);

void run_server(long port, long n1, long n2, long n3)
{
    using connect_t = std::unique_ptr<tcp::iostream>;
    Queue<connect_t> queue;
    
    auto prgm = [&queue](long n1, long n2, long n3) {
        for (;;) {
            connect_t conn = std::move(queue.pop());
            if (!conn)
                continue;
            play_server(*conn, n1, n2, n3);
            conn->close();
        }
    };

    std::vector<std::thread> workers;
    for (int w = 0; w < 30; ++w) { // 30 workers
        workers.emplace_back(prgm, n1, n2, n3);
    }

    boost::asio::io_service ios;
    tcp::endpoint endpoint(tcp::v4(), port);
    tcp::acceptor acceptor(ios, endpoint);
    for (;;) {
        boost::system::error_code err;
        tcp::iostream *conn= new tcp::iostream();
        acceptor.accept(*(conn->rdbuf()), err);
        if (!err) {
            connect_t _conn(conn);
            queue.push(std::move(_conn));
        } else {
            delete conn;
        }
    }
}

int run_client(std::string const& addr, long port, long n1, long n2, long n3) {
	dgk_prvkey_t *sk;
	dgk_pubkey_t *pk;
	dgk_keygen(1024, 16, &pk, &sk); // 16-bit plaintexts
    tcp::iostream conn(addr, std::to_string(port));
    if (!conn) {
        std::cerr << "Can not connect to server!" << std::endl;
        return -1;
    }
    NetworkLog::bytes_sent = 0;
    NetworkLog::bytes_recev = 0;
    double all_time = 0.;
    do {
        AutoTimer time(&all_time);
        play_client(conn, sk, pk, n1, n2, n3);
    } while(0);
    printf("%.3f\n", all_time);
    conn.close();
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

    EncMatrix enc_A(n1 * n2);
	for (long r = 0; r < n1; r++) {
		long offset = r * n2;
		for (long c = 0; c < n2; c++) {
			mpz_init(enc_A.at(offset + c));
			receive_mpz(enc_A[offset + c], buffer, conn);
		}
	}
	
	mpz_t tmp, zero;
	mpz_init(tmp);
	mpz_init(zero);
	mpz_set_ui(zero, 0);

	for (long r = 0; r < n1; r++) {
		long offset = r * n2;
		for (long c = 0; c < n3; c++) {
			mpz_t res;
			mpz_init(res);
			dgk_encrypt_plain(res, &pk, zero, gmp_rand);
			for (long k = 0; k < n2; k++) {
				dgk_hom_mult(tmp, enc_A.at(offset + k), B[k][c], &pk);
				dgk_hom_add(res, res, tmp, &pk);
			}
			send_mpz(res, buffer, conn);
		}
	}
	mpz_clear(tmp);
}

void play_client(std::iostream &conn, dgk_prvkey_t *sk, dgk_pubkey_t *pk, long n1, long n2, long n3)
{
	Matrix A, B;
	NTL::SetSeed(NTL::to_ZZ(123)); // for debugging
	A.SetDims(n1, n2);
	B.SetDims(n2, n3);
	randomize_matrix(A);
	randomize_matrix(B);

	std::vector<uint32_t> buffer(2048 >> 5);
	send_pk(pk, buffer, conn);
	
	EncMatrix enc_A(n1 * n2);
	encrypt_matrix(enc_A, A, pk);
	for (long r = 0; r < n1; r++) {
		int offset = r * n2;
		for (long c = 0; c < n2; c++) {
			send_mpz(enc_A.at(offset + c), buffer, conn);
		}
	}
	
	mpz_t plt, ctx;
	mpz_init(plt);
	mpz_init(ctx);
	for (long r = 0; r < n1; r++) {
		for (long c = 0; c < n3; c++) {
			receive_mpz(ctx, buffer, conn);
			dgk_decrypt(plt, pk, sk, ctx);
		}
	}
	mpz_clear(plt);
	mpz_clear(ctx);
}
