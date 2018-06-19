#include <boost/asio.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <thread>
#include <iostream>
#include <chrono>
using boost::asio::ip::tcp;

int run_server(long port) {
    boost::asio::io_service ios;
    tcp::endpoint endpoint(tcp::v4(), port);
    tcp::acceptor acceptor(ios, endpoint);
    std::atomic<int> id(0);
    auto serve = [&id, &acceptor]() {
        for (;;) {
            tcp::iostream conn;
            boost::system::error_code err;
            acceptor.accept(*conn.rdbuf(), err);

            if (!err) {
                int nxt = id.fetch_add(1);
                conn << "SERVING" << std::to_string(nxt);
                conn.close();
            }
        }
    };

    std::vector<std::thread> workers;
    for (int t = 0; t < 4; ++t)
        workers.emplace_back(serve);

    for (auto && w: workers)
        w.join();
    return 0;
}

int run_client(std::string addr, long port) {
    tcp::iostream conn(addr, std::to_string(port));
    if (conn) {
        std::string resp;
        conn >> resp;
        conn.close();
        std::cout << resp << '\n';
    } else {
        std::cerr << "Can not connect:" << '\n';
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        run_server(12346);
    } else {
        int rep = std::atoi(argv[1]);
        std::cout << "to repeat " << rep << '\n';
        std::vector<std::thread> workers;
        std::string addr{"127.0.0.1"};
        for (int i = 0; i < rep; ++i)
            workers.emplace_back(run_client, addr, 12346);

        for (auto &&w : workers)
            w.join();
    }
    return 0;
}
