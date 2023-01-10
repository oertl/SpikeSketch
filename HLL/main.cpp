#include <iostream>
#include <random>

#include "HyperLogLog.hpp"
#include <fstream>


using namespace std;

int main() {
    int registerNum = 108;
    int nf = 1 << 10;//Actual Cardinality

    srand(time(NULL));
    uint64_t key;
    int memory =registerNum*6;
    time_t seed = rand() + 1316561;
    mt19937 rng(seed);
    uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    hyperlogloglog::HyperLogLog<> hll(registerNum);
    for (key = 1; key <= nf; key++) {
        uint64_t y = dist(rng);
        hll.add(y);
    }
    double result = hll.estimate();

    char traceFilePath[100];
    ofstream logFile;
    sprintf(traceFilePath, "result.txt");
    logFile.open(traceFilePath, ios::out);
    logFile <<"Input:\n\tReal Cardinality: "<<nf<< "\n\tNumber of Registers: " << registerNum << endl;
    logFile << "######################################" << endl;
    logFile <<"Output:\n\t"<<"Estimate Cardinality: "<<result<< "\n\tMemory Usage: " << memory;
    logFile.close();

    return 0;
}
