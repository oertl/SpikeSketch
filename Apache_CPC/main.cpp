#include <iostream>
#include "cpc_sketch.hpp"
#include <fstream>

using namespace std;

int main() {
    uint64_t nf = 1 << 25;
    int regNum = 49;

    srand(time(NULL));
    uint64_t key;
    double value = 0;

    time_t seed = rand() + 1316561;
    mt19937 rng(seed);
    uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    datasketches::cpc_sketch_alloc<std::allocator<uint8_t>> cpc_sk(regNum);
    double maxSize = 0.0;
    for (key = 1; key <= nf; key++) {
        uint64_t y = dist(rng);
        uint8_t lg_size = cpc_sk.update(y);
        double size = (double) (1 << lg_size);
        if (size > maxSize) {
            maxSize = size;
        }
        double result = cpc_sk.get_estimate();
        value = result;
    }
    double Mem = 8 * regNum + (ceil(std::log2(regNum)) + 6) * maxSize + 6;

    ofstream logFile;
    logFile.open("result.txt", ios::out);
    logFile << "Input:\n\t" << "Real Cardinality: " << nf << "\n\tNumber of Register: " << regNum << endl;
    logFile << "######################################" << endl;
    logFile << "Output:\n\t" << "Estimate Cardinality: " << value << "\n\tMemory Usage: " << Mem;
    logFile.close();
    return 0;
}
