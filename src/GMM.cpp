/// General Matrix Multiplication
#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/NumbTh.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Matrix.hpp"
#include "CryptGMM/Timer.hpp"
#include "CryptGMM/HElib.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <iostream>
#include <numeric>
#include <list>
using boost::asio::ip::tcp;

std::pair<double, double> mean_std(std::vector<double> const& times) {
    if (times.empty())
        return {0., 0.};
    if (times.size() == 1)
        return {times[0], 0.};
    double mean = 0.;
    for (double t : times)
        mean += t;
    mean /= times.size();
    double stdr = 0.;
    for (double t : times)
        stdr += (t - mean) * (t - mean);
    stdr = std::sqrt(stdr / (times.size() - 1));
    return {mean, stdr};
}

inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

void zero(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = 0;
}

void randomize(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = j;
    //NTL::RandomBnd(10L);//i * 10 + j;
}

struct Duplication {
    long blk_id;
    long dup;
    long slots_left;
};

template <class T>
void rotate_vector(T &vec, long k) {
    long n = vec.size();
    T tmp(vec);
    for (long i = 0; i < n; i++)
        tmp[((i+k)%n + n)%n] = vec[i];
    vec = tmp;
}

Duplication compute_duplications(const long total_rows,
                                 const long x_blk_id,
                                 const long num_slots) {
    long N = std::min((x_blk_id  + 1) * num_slots, total_rows) - x_blk_id * num_slots;
    Duplication d;
    d.blk_id = x_blk_id;
    d.dup = num_slots / N;
    d.slots_left = num_slots - d.dup * N;
    return d;
}

std::vector<long> compute_rotations(long x_blk, long total_rows_X,
                                    long y_blk, long total_rows_Y,
                                    long num_slots) {
    auto dupA = compute_duplications(total_rows_X, x_blk, num_slots);
    auto dupB = compute_duplications(total_rows_Y, y_blk, num_slots);
    std::vector<long> rotations;
    long sze;
    // fully packed
    if (dupA.dup == 1 and dupB.dup == 1) {
        sze = num_slots;
    } else if (dupA.slots_left != 0 or dupB.slots_left != 0) {
        // cumble to implement, so just skip this case
        std::cout << "cumble\n" << std::endl;
        sze = num_slots;
    } else {
        long NA = (num_slots - dupA.slots_left) / dupA.dup;
        long NB = (num_slots - dupB.slots_left) / dupB.dup;
        sze = std::max(NA, NB) / (dupA.dup * dupB.dup);
        sze = std::max(1L, sze);
    }
    rotations.resize(sze);
    std::iota(rotations.begin(), rotations.end(), 0);
    return rotations;
}

void fill_compute(Matrix& mat,
                  int x, int y, int k,
                  const std::vector<NTL::zz_pX> &slots,
                  const EncryptedArray *ea) {
    const long l = ea->size();
    const long d = ea->getDegree();
    long row_start = x * l;
    long row_end = row_start + l;
    long col_start = y * l;
    long col_end = col_start + l;

    assert(slots.size() == l);
    for (long ll = 0; ll < l; ll++) {
        long computed = NTL::coeff(slots[ll], d - 1)._zz_p__rep;
        long row = row_start + ll;
        long col = col_start + ll + k;
        if (col >= col_end)
            col = col_start + col % l;
		col %= mat.NumCols();
		row %= mat.NumRows();
		mat.put(row, col, computed);
    }
}

long ceil_round(long a, long l) {
    return (a + l - 1) / l * l;
}

FHEcontext receive_context(std::istream &s) {
    unsigned long m, p, r;
    std::vector<long> gens, ords;
    readContextBase(s, m, p, r, gens, ords);
    FHEcontext context(m, p, r, gens, ords);
    NTL::zz_p::init(p);
    s >> context;
    return context;
}

void send_context(std::ostream &s, FHEcontext const& context) {
    writeContextBase(s, context);
    s << context;
}

struct ClientBenchmark {
    std::vector<double> enc_times;
    std::vector<double> dec_times;
    std::vector<double> total_times;
};
ClientBenchmark clt_ben;

void play_client(tcp::iostream &conn,
                 FHESecKey &sk,
                 FHEcontext &context,
                 const long n1,
                 const long n2,
                 const long n3,
                 bool verbose) {
    FHEPubKey ek(sk);
    ek.makeSymmetric();
    conn << ek;
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();

    NTL::SetSeed(NTL::to_ZZ(123));
    Matrix A, B, ground_truth;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A);
    randomize(B);
    ground_truth = mul(A, B);
    /// print grouth truth for debugging
    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
    const long MAX_X2 = round_div(B.NumCols(), l);

    std::vector<std::vector<Ctxt>> uploading;
    uploading.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, sk));
	double enc_time;
	do {
		AutoTimer timer(&enc_time);
    	/// encrypt matrix
		NTL::ZZX packed_poly;
		for (int x = 0; x < MAX_X1; x++) {
			for (int k = 0; k < MAX_Y1; k++) {
				internal::BlockId blk = {x, k};
				auto block = internal::partition(A, blk, *ea, false);
				rawEncode(packed_poly, block.polys, context);
				sk.Encrypt(uploading[x][k], packed_poly);
			}
		}
	} while(0);
    clt_ben.enc_times.push_back(enc_time);

	double sent_ctx_time;
	do {
		AutoTimer timer(&sent_ctx_time);
		/// send ciphertexts of matrix
		for (auto const& row : uploading) {
			for (auto const& ctx : row)
				conn << ctx;
		}
	} while (0);
    if (verbose)
        std::cout << "Sent " << MAX_X1 * MAX_Y1 << " ctxts in " << sent_ctx_time << " ms" << std::endl;

    /// waiting results
    long rows_of_A = A.NumRows();
    long rows_of_Bt = B.NumCols(); // Bt::Rows = B::Cols
    std::list<Ctxt> ret_ctxs;
    for (int x = 0; x < MAX_X1; x++) {
        for (int y = 0; y < MAX_X2; y++) {
            size_t num_rots = compute_rotations(x, rows_of_A,
                                                y, rows_of_Bt, l).size();
            for (size_t k = 0; k < num_rots; k++) {
                Ctxt result(sk);
                conn >> result;
                ret_ctxs.emplace_back(result);
            }
        }
    }
    if (verbose)
        std::cout << "Received " << ret_ctxs.size() << " ctxts" << std::endl;
    /// decrypt
    Matrix computed;
    computed.SetDims(A.NumRows(), B.NumCols());
    zero(computed);
    int x = 0;
    int y = 0;
    auto itr = ret_ctxs.begin();
    std::vector<NTL::zz_pX> slots;
    NTL::ZZX decrypted;
	double decrypt_time;
	do {
		AutoTimer timer(&decrypt_time);
		for (int x = 0; x < MAX_X1; x++) {
			for (int y = 0; y < MAX_X2; y++) {
				size_t num_rots  = compute_rotations(x, rows_of_A,
													 y, rows_of_Bt, l).size();
				for (size_t k = 0; k < num_rots; k++) {
					sk.Decrypt(decrypted, *itr++);
					rawDecode(slots, decrypted, context);
					fill_compute(computed, x, y, k, slots, ea);
				}
			}
		}
	} while (0);
    clt_ben.dec_times.push_back(decrypt_time);
	if (!::is_same(ground_truth, computed, NTL::zz_p::modulus()))
	  	std::cerr << "The computation seems wrong " << std::endl;
}

struct ServerBenchmark {
    std::vector<double> eval_times;
};
ServerBenchmark srv_ben;

void play_server(tcp::iostream &conn,
                 const long n1,
                 const long n2,
                 const long n3) {
    auto context = receive_context(conn);
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    FHEPubKey ek(context);
    conn >> ek;
    NTL::SetSeed(NTL::to_ZZ(123));

    Matrix A, B;
    A.SetDims(n1, n2);
    B.SetDims(n2, n3);
    randomize(A);
    randomize(B);
    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
    Matrix Bt;
    transpose(&Bt, B);
    const long MAX_X2 = round_div(Bt.NumRows(), l);
    const long MAX_Y2 = round_div(Bt.NumCols(), d);
    NTL::Mat<internal::PackedRows> precomputed; // 2D-array of polynomials
    precomputed.SetDims(MAX_X2, MAX_Y1);
    for (int y = 0; y < MAX_X2; y++) {
        for (int k = 0; k < MAX_Y1; k++) {
            internal::BlockId blk = {y, k};
            precomputed[y][k] = internal::partition(Bt, blk, *ea, true);
        }
    }

    /// receving ciphertexts from the client
    std::vector<std::vector<Ctxt>> received;
    received.resize(MAX_X1, std::vector<Ctxt>(MAX_Y1, ek));
    for (int x = 0; x < MAX_X1; x++) {
        for (int k = 0; k < MAX_Y1; k++)
            conn >> received[x][k];
    }
    /// compute the matrix mulitplication
    long rows_of_A = A.NumRows();
    long rows_of_Bt = Bt.NumRows();
    size_t send_ctx = 0;
    NTL::ZZX packed_polys;
    double computation{0.}, network{0.};
	for (int x = 0; x < MAX_X1; x++) {
		for (int y = 0; y < MAX_X2; y++) {
			std::vector<long> rotations = compute_rotations(x, rows_of_A,
															y, rows_of_Bt, l);
			size_t num_rotations = rotations.size();
			std::vector<Ctxt> summations(num_rotations, ek);
			double eval_time;
			do {
				AutoTimer timer(&eval_time);
				for (int k = 0; k < MAX_Y1; k++) {
					for (size_t rot = 0; rot < num_rotations; rot++) {
						auto rotated(precomputed[y][k].polys);
						rotate_vector(rotated, -rotations[rot]);
						/// pack the polys into one poly
						rawEncode(packed_polys, rotated, context);
						Ctxt client_ctx(received[x][k]);
						/// ctx * plain
						client_ctx.multByConstant(packed_polys);
						summations[rot] += client_ctx;
					}
				}
				for (auto &sm : summations)
					sm.modDownToLevel(1);
			} while (0);
			computation += eval_time;

			// start_time = Clock::now();
			double send_ctx_time;
			do {
				AutoTimer timer(&send_ctx_time);
				for (auto const&sm : summations)
					conn << sm;
			} while (0);
			network += send_ctx_time;
			send_ctx += summations.size();
		}
	}
    srv_ben.eval_times.push_back(computation);
}

int run_client(long n1, long n2, long n3, bool verbose) {
    tcp::iostream conn("192.168.0.116", "12345");
    if (!conn) {
        std::cerr << "Can not connect to server!" << std::endl;
        return -1;
    }
    const long m = 8192;
    //const long p = 401;
    const long p = 769;
    const long r = 2;
    const long L = 4;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 15;// + std::ceil(std::log(m)/2 + r * std::log(p));
    buildModChain(context, L);
    if (verbose) {
        std::cout << "kappa = " << context.securityLevel() << std::endl;
        std::cout << "slot = " << context.ea->size() << std::endl;
        std::cout << "degree = " << context.ea->getDegree() << std::endl;
    }
    FHESecKey sk(context);
    sk.GenSecKey(64);
    /// send FHEcontext obj
    send_context(conn, context);
    double all_time;
    do {
        AutoTimer time(&all_time);
        /// send the evaluation key
        play_client(conn, sk, context, n1, n2, n3, verbose);
    } while(0);
    clt_ben.total_times.push_back(all_time);
    conn.close();
    return 1;
}

int run_server(long n1, long n2, long n3) {
    boost::asio::io_service ios;
    tcp::endpoint endpoint(tcp::v4(), 12345);
    tcp::acceptor acceptor(ios, endpoint);

    for (long run = 0; run < 20; run++) {
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
    argmap.arg("N", n1, "n1");
    argmap.arg("M", n2, "n2");
    argmap.arg("D", n3, "n3");
    argmap.arg("R", role, "role. 0 for server and 1 for client");
    argmap.parse(argc, argv);
    if (role == 0) {
        run_server(n1, n2, n3);
        auto eval_time = mean_std(srv_ben.eval_times);
        printf("Server evaluation time\n");
        printf("%.3f \\pm %.3f\n", eval_time.first, eval_time.second);
    } else if (role == 1) {
        for (long run = 0; run < 20; run++)
            run_client(n1, n2, n3, run == 0);
        printf("Client enc dec total\n");
        auto time = mean_std(clt_ben.enc_times);
        printf("%.3f \\pm %.3f ", time.first, time.second);
        time = mean_std(clt_ben.dec_times);
        printf("%.3f \\pm %.3f ", time.first, time.second);
        time = mean_std(clt_ben.total_times);
        printf("%.3f \\pm %.3f\n", time.first, time.second);
    } else {
		argmap.usage("General Matrix Multiplication for |N*M| * |M*D|");
		return -1;
	}
}
