/// General Matrix Multiplication
#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/NumbTh.h>

#include "SMP/Matrix.hpp"
#include "SMP/Timer.hpp"
#include "SMP/HElib.hpp"
#include "SMP/literal.hpp"
#include "SMP/network/net_io.hpp"
#include "SMP/SMPServer.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <queue>
#include <iostream>
#include <thread>
#include <future>
#include <numeric>
#include <list>
#include <memory>
using boost::asio::ip::tcp;
constexpr int REPEAT = 1;

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
            mat[i][j] = NTL::RandomBnd(p);
}


void fill_compute(Matrix& mat,
				  long row_blk,
				  long col,
                  const std::vector<long> &inner_prod,
                  const EncryptedArray *ea)
{
    const long l = ea->size();
    assert(inner_prod.size() == l);
	const bool is_vec = mat.NumRows() == 1;
	const long row_start = is_vec ? 0 : row_blk * l;
	const long col_start = is_vec ? row_blk * l : col;

    for (long ll = 0; ll < l; ll++) {
        long computed = inner_prod[ll];
		if (!is_vec) {
			long row = row_start + ll;
			if (row < mat.NumRows())
				mat.put(row, col, computed);
			else
				break;
		} else {
			long col = col_start + ll;
			if (col < mat.NumCols())
				mat.put(0, col, computed);
			else
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

    Matrix A, B, ground_truth;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    NTL::SetSeed(NTL::to_ZZ(123));
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
			double one_pack_time = 0.; 
			double one_enc_time = 0.;
			auto block = internal::partition(A, blk, ea, false);
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
	conn.flush();
    clt_ben.ctx_sent = MAX_X1 * MAX_Y1;
	/// we convert DoubleCRT to poly form when send ciphertexts through, and thus, we
	/// count this cost as a part of encryption.
	clt_ben.enc_times.back() += getTimerByName("TO_POLY_OUTPUT")->getTime() * 1000.;

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
    //NTL::Vec<long> decrypted;
	NTL::ZZX decrypted;
	double decrypt_time = 0.;
    double unpack_time = 0.;
	long ctx_idx = 0;
	bool dec_pass = true;
	for (const auto &ctx : ret_ctxs) {
		double one_dec_time = 0.; 
		double one_unpack_time = 0.;
		do {
			AutoTimer timer(&one_dec_time);
			dec_pass &= ctx.isCorrect();
			//faster_decrypt(decrypted, sk, ctx);
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
	/// we convert poly to DoubleCRT when receiving ciphertexts.
	decrypt_time += getTimerByName("FROM_POLY_OUTPUT")->getTime() * 1000.;
    clt_ben.dec_times.push_back(decrypt_time);
    clt_ben.unpack_times.push_back(unpack_time);
	if (!::is_same(ground_truth, computed, NTL::zz_p::modulus()))
		std::cerr << "The computation seems wrong " << std::endl;
	if (!dec_pass)
		std::cerr << "Decryption might fail" << std::endl;
}

int run_client(std::string const& addr, long port,
               long n1, long n2, long n3) {
    const long m = 8192;
    const long p = 84961;//70913;
    const long r = 1;
    const long L = 2;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 60;
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
        double all_time = 0.;
        do {
            send_context(conn, context);
            AutoTimer time(&all_time);
            /// send the evaluation key
            play_client(conn, sk, context, n1, n2, n3);
        } while(0);
        clt_ben.total_times.push_back(all_time);
        conn.close();
		resetAllTimers(); // reset timers in HElib
    }
    return 1;
}

template <typename T>
class Queue
{
public:

    T pop()
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        while (queue_.empty())
        {
            cond_.wait(mlock);
        }
        auto item = std::move(queue_.front());
        queue_.pop();
        return item;
    }

    void pop(T& item)
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        while (queue_.empty())
        {
            cond_.wait(mlock);
        }
        item = std::move(queue_.front());
        queue_.pop();
    }

    void push(const T& item)
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        queue_.push(item);
        mlock.unlock();
        cond_.notify_one();
    }

    void push(T&& item)
    {
        std::unique_lock<std::mutex> mlock(mutex_);
        queue_.push(std::move(item));
        mlock.unlock();
        cond_.notify_one();
    }

private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_;
};

int run_server(long port, long n1, long n2, long n3) {
    using connect_t = std::unique_ptr<tcp::iostream>;
    Queue<connect_t> queue;
    
    auto prgm = [&queue](long n1, long n2, long n3) {
        for (;;) {
            connect_t conn = std::move(queue.pop());
            if (!conn)
                continue;
            SMPServer server;
            server.run(*conn, n1, n2, n3);
            conn->close();
        }
    };

    std::vector<std::thread> workers;
    for (int w = 0; w < 60; ++w) {
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
            //std::cout << "Receive from " << conn->rdbuf()->remote_endpoint().address().to_string() << "@";
            //std::cout << conn->rdbuf()->remote_endpoint().port() << "\n";
            connect_t _conn(conn);
            queue.push(std::move(_conn));
        } else {
            delete conn;
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
