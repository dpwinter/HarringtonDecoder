#ifndef TORICCODE_H_
#define TORICCODE_H_

#include <random>

class ToricCode {
    private:
        int L;
        bool*** qubits;
        bool** stabs;

        std::random_device randDev; // wraps /dev/urandom
        std::mt19937 randGen{randDev()}; // init w/ random seed
        std::uniform_real_distribution<double> randDist;

    public:
        ToricCode(int L);
        virtual ~ToricCode();
        void reset();
        void flip(int i, int j, int dir);
        bool getStab(int i, int j);
        bool** getSyndromes();
        bool getQubit(int i, int j, int k);
        bool hasLogErr();
        void noise(double p);
        void setSeed(int seed) { this->randGen.seed(seed); };
};

#endif
