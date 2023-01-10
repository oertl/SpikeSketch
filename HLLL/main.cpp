#include <iostream>
#include <random>
#include "HyperLogLogLog.hpp"
#include <fstream>

using namespace std;

int main() {
    uint64_t nf = 1 << 10;
    int regNum = 109;

    srand(time(NULL));
    uint64_t key;
    time_t seed = rand() + 1316561;
    mt19937 rng(seed);
    uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    hyperlogloglog::HyperLogLogLog<> hlll(regNum);
    double maxSize = 0.0;
    for (key = 1; key <= nf; key++) {
        uint64_t y = dist(rng);
        double size = hlll.add(y);
        if (size > maxSize) {
            maxSize = size;
        }
    }
    double result = hlll.estimate();
    double Mem = 3.0 * regNum + (ceil(std::log2(regNum)) + 6) * maxSize + 6;

    ofstream logFile;
    logFile.open("result.txt", ios::out);
    logFile << "Input:\n\t" << "Actual Cardinality: " << nf << "\n\tNumber of Register: " << regNum << endl;
    logFile << "######################################" << endl;
    logFile << "Output:\n\t" << "Estimate Cardinality: " << result << "\n\tMemory Usage: " << Mem;
    logFile.close();
    return 0;































//
//    for (int nfi=0;nfi<numOftestCards;nfi++)
//    {
//        uint64_t nf = nfs[nfi];
////        logFile<<"processing "<<nfi<<"cards: "<<nf<<endl;
//
//        double avgSize=0.0;
//        double totalMaxSize=0.0;
//        for (int roundIdx = 0; roundIdx < round; roundIdx++){
//
//            time_t seed = rand()+1316561;
//            mt19937 rng(seed);
//            uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
//            hyperlogloglog::HyperLogLogLog<> hlll (registerNum);
//
//
//            double maxSize=0.0;
//            for (key = 1; key <= nf; key++) {
//                uint64_t y = dist(rng);
//                double size=hlll.add(y);
//                if(size>maxSize){
//                    maxSize=size;
//                }
//            }
//
//            avgSize+= maxSize;
//            if(maxSize>totalMaxSize){
//                totalMaxSize=maxSize;
//            }
//
//            double result= hlll.estimate();
//            values[nfi][roundIdx]= result;
//
////            if((roundIdx+1)%(1)==0){
////                logFile<<nf<<": "<<(roundIdx+1)<<"rounds, avgSize: "<<avgSize/(roundIdx+1.0)<<", maxSize: "<<totalMaxSize<<endl;
////            }
//        }
//        avgSize=avgSize/ round;
//
//
//        double mean = stats_mean(values[nfi], round);
//        double coe=nf/mean;
//
//        double totalAvgMem=3*registerNum+(hyperlogloglog::log2i(registerNum)+6)*avgSize;
//        double totalMaxMem=3*registerNum+(hyperlogloglog::log2i(registerNum)+6)*totalMaxSize;
//
//        logFile<<"nf: "<<nf<<", mean: "<<mean<<", coe: "<<coe<<", avgSize: "<<avgSize<<", maxSize: "<<totalMaxSize<<", totalAvgMem: "<<totalAvgMem<<", totalMaxMem: "<<totalMaxMem<<endl;//实际基数
//        coeSum+=coe;
//
//
//        for (int roundIdx = 0; roundIdx < round; roundIdx++){
//            values[nfi][roundIdx] =values[nfi][roundIdx]/nf;
//        }
//        double varer =stats_variance(values[nfi],1.0, round);
//        double stdErr=sqrt(varer);
//        double MVP=varer*memoryBits;
//        logFile<<"nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
//        dataFile<<stdErr<<","<<flush;
//    }
//    dataFile<<"],"<<endl;
//    dataFile<<"]"<<endl;
//
//    logFile.close();
//    dataFile.close();
//
//    return 0;


//    mt19937 rng(0);
//    uniform_int_distribution<uint64_t> dist(0,0x3ff);
//    hyperlogloglog::HyperLogLog<> hll(16);
//    for (int i = 0; i < 1000; ++i) {
//        uint64_t y = dist(rng);
//        hll.add(y);
//    }
//    double result=hll.estimate();
//    cout << result << endl;
//    return 0;
}
