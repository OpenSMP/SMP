//
// Created by riku on 1/3/18.
//
#include <NTL/ZZ.h>
#include <HElib/NumbTh.h>
#include <iostream>
void MiniONN(long m) {
    long p = NTL::RandomPrime_long(10, 20);
    while (true) {
        long d = multOrd(p, m);
        if (d == 1)
            break;
        do {
            p = NTL::NextPrime(p + 2);
        } while (m % p == 0);
    };
    std::cout << m << " " << p << std::endl;
}

int main() {
    MiniONN(8192);
    return 0;
}

