#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <iostream>
#include <fstream>

void run_client() {
    FHEcontext context(1024, 1031, 1);
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);
    std::ofstream of("pk.pk", std::ios::binary);
    FHEPubKey pk(sk);
    // pk.makeSymmetric();
    of << pk;
    of.close();

    Ctxt ctx(sk);
    sk.Encrypt(ctx, NTL::to_ZZX(3));

    using boost::asio::ip::tcp;
    tcp::iostream connect("127.0.0.1", "12345");
    if (!connect) {
        std::cerr << "Can not connect:" << 
             connect.error().message() << std::endl;
        return;
    }
    writeContextBase(connect, context);
    connect << context;
    connect << ctx;
    connect.flush();
    connect >> ctx;
    NTL::ZZX dec;
    sk.Decrypt(dec, ctx);
    std::cout << "return " << dec[0] << std::endl;
    connect.close();
}

void run_server() {
    using boost::asio::ip::tcp;
    boost::asio::io_service ios;
    tcp::endpoint endpoint(tcp::v4(), 12345);
    tcp::acceptor acceptor(ios, endpoint);

    for (;;) {
        tcp::iostream client;
        boost::system::error_code err;
        acceptor.accept(*client.rdbuf(), err);
        if (!err) {
            std::cout << "connected" << std::endl;
            unsigned long m, p, r;
            std::vector<long> gens, ords;
            readContextBase(client, m, p, r, gens, ords);
            FHEcontext context(m, p, r, gens, ords);
            client >> context;
            std::ifstream ifs("pk.pk", std::ios::binary);
            FHEPubKey pk(context);
            ifs >> pk;
            ifs.close();

            Ctxt ctx(pk);
            client >> ctx;
            ctx.multiplyBy(ctx);
            ctx.modDownToLevel(1);
            client << ctx;
            client.close();
            break;
        } else {
            std::cerr << "Error " << err << std::endl;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        run_client();
    } else {
        run_server();
    }
}
