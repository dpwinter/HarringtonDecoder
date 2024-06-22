#include "Location.h"
#include "ToricCode.h"
#include "CA.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <fstream>

class Timer {
    using clk = std::chrono::steady_clock;
    clk::time_point begin;
    clk::time_point end;

    public:
        void clear() { begin = end = clk::now(); }
        void start() { begin = clk::now(); }
        void stop() { end = clk::now(); }

        friend std::ostream& operator<<(std::ostream& o, const Timer& timer) {
            return o << timer.secs(); // string repr
        }

        double secs() const {
            if(end <= begin) { return 0.0; }
            auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
            return diff.count() / 1000000.0;
        }
};

double mean(std::vector<int> &counts, int N) {
    double sum = 0;
    for(int c : counts) {
        sum += c;
    }
    return sum / N;
}

double sd(std::vector<int> &counts, int N) {
    int n = counts.size();
    double mu = mean(counts, N);
    double sqDistSum = 0;
    for(int c : counts) {
        double dist = mu - c;
        sqDistSum += std::pow(dist, 2);
    }
    return std::sqrt( sqDistSum / N );
}

void applyCorrections(ToricCode& tc, Location** corrections, int L) {
    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++) {
            tc.flip(i,j,corrections[i][j]);
        }
    }
}

double harringtonVis(ToricCode &tc, CA &ca, double p, int N, int L) {
    std::ofstream qubits_file;
    std::ofstream flips_file;
    std::ofstream counts_file;
    qubits_file.open("qubits.csv", std::ios_base::out); // overwrite
    flips_file.open("flipsigs.csv", std::ios_base::out); // overwrite
    counts_file.open("countsigs.csv", std::ios_base::out); // overwrite

    tc.reset();
    ca.reset();

    for(int i=0; i<N; i++) {
        tc.noise(p);

        // qubits
        for(int j=0; j<L; j++){
            for(int k=0; k<L; k++) {
                for (int l=0; l<2; l++) {
                    if (tc.getQubit(j,k,l)) {
                        qubits_file << j << "," << k << "," << l << " ";
                    }
                }
            }
        }
        qubits_file << std::endl;

        bool** syndromes = tc.getSyndromes();
        Location** corrections = ca.step(syndromes);

        // signals
        for(int j=0; j<L; j++){
            for(int k=0; k<L; k++) {
                int n_counts = 0;
                for(int l=0; l<8; l++) {
                    if(ca.getCell(j,k)->getMemory(0)->flipSig[l]) { // TODO: include flipSigs on all levels.
                        flips_file << j << "," << k << "," << l << " ";
                    }
                    n_counts += ca.getCell(j,k)->getMemory(0)->countSig[l];
                }
                if(n_counts > 0) {
                    counts_file << j << "," << k << "," << 0 << " ";
                }
            }
        }

        flips_file << std::endl;
        counts_file << std::endl;

        applyCorrections(tc, corrections, L);

    }
    return 0.0;
}

double benchmarkHarrington(ToricCode &tc, CA &ca, double p, int N, int L) {

    std::ofstream counts_file;
    counts_file.open("./data/L=" + std::to_string(L) + "_p=" + std::to_string(int(1/p)) + ".csv", std::ios_base::app); // append

    int tot_count = 0;

    for(int n=0; n<N; n++) {
        tc.reset();
        ca.reset();

        int count = 0;

        while(!tc.hasLogErr()) {
            tc.noise(p);
            bool** syndromes = tc.getSyndromes();
            Location** corrections = ca.step(syndromes);

            applyCorrections(tc, corrections, L);

            count += 1;
        }
        counts_file << std::to_string(count) << std::endl;
        tot_count += count;
    }
    return tot_count;
}

// main loop //

int main() {

    // Q=3 (hard-coded)
    int U = 10;
    double fC = 9/10.;
    double fN = 4/10.;
    int N = 100;

    std::vector<int> Ls = {9};
    // std::vector<double> ps = {1e-1,5e-2,1e-2,5e-3,3e-3,2e-3};
    // std::vector<double> ps = {11,12,14,16,25,33,50,111,125,142,166,250};
    std::vector<double> ps = {3e-3};



    Timer timer;
    timer.start();

    for(int L : Ls) {

        std::cout << "--- Lattice size " << L << " ---\n";
        std::vector<int> counts(ps.size(), 0);

        ToricCode tc(L);
        CA ca(L,U,fC,fN);

        for(int i=0; i<ps.size(); i++) {
            counts[i] = benchmarkHarrington(tc, ca, ps[i], N, L);
            // counts[i] = benchmarkToricCode(tc, ps[i], N);
            // counts[i] = harringtonVis(tc, ca, ps[i], N, L);
        }

        for(int i=0; i<ps.size(); i++) {
            std::cout << "p=" << ps[i] << ": mu=" << static_cast<double>(counts[i]) / N << '\n';
        }
    }

    timer.stop();
    std::cout << "[[ time: " << timer << " secs ]]\n";

    return 0;
}
