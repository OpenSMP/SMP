#ifndef CRYPT_GMM_SMP_SERVER_HPP
#define CRYPT_GMM_SMP_SERVER_HPP
#include <HElib/EncryptedArray.h>
#include <HElib/Ctxt.h>
#include <HElib/NumbTh.h>

#include "SMP/Matrix.hpp"
#include "SMP/DoublePacking.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <vector>
#include <atomic>
#include <thread>
#include <list>
class FHEcontext;
using boost::asio::ip::tcp;
class SMPServer {
public:
	SMPServer();

	~SMPServer();

	void run(tcp::iostream &conn,
			 const long n1,
			 const long n2,
			 const long n3);

protected:
	void setup(tcp::iostream &conn);

	void process_columns();

	void receive_ctx(tcp::iostream &conn);

	void evaluate();

	void evaluate_mat_vec();

	void extract_and_merge();

	void response_ctx(tcp::iostream &conn);

private:
	Matrix A, B;
	Matrix ground_truth;

	FHEcontext *context = nullptr;
	FHEPubKey *ek = nullptr;
	std::vector<std::vector<Ctxt>> enc_A_blk;
	NTL::Mat<internal::PackedRows> plain_B_blk;
	NTL::Mat<zzX> encoded_plain_B_blk;
	std::vector<std::vector<zzX>> linear_map_coeffs;
	std::list<Ctxt> results;
};
#endif // CRYPT_GMM_SMP_SERVER_HPP
