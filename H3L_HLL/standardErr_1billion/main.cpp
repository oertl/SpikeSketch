#include <iostream>
#include <random>

#include "HyperLogLog.hpp"
#include <random>
#include <fstream>
#include "stats.h"
using namespace std;

int main() {
    srand(time(NULL));

    int round = 100;

    vector<int> memorys={648,2592,16200};

    int numOfMemValues=memorys.size();
    double** values=new double*[numOfMemValues];
    for(int i=0;i<numOfMemValues;i++){
        values[i]=new double[round];
    }

    uint64_t key;
    int nf=1<<30;

    char traceFilePath[100];
    ofstream logFile;
    sprintf(traceFilePath, "total.log");
    logFile.open(traceFilePath,ios::out);

    ofstream dataFile;
    sprintf(traceFilePath, "total.dat");
    dataFile.open(traceFilePath,ios::out);


    logFile<<"nf: "<<nf<<endl;//实验次数
    logFile<<"round: "<<round<<endl;//实验次数
    logFile<<"######################################\n\n\n"<<endl;


    dataFile<<"[";
    dataFile<<endl;

    dataFile<<"[";
    for (int memIdx=0;memIdx<numOfMemValues;memIdx++)
    {
        int bits = memorys[memIdx];
        dataFile<<bits<<",";
    }
    dataFile<<"],"<<endl;

    double coeSum=0.0;


    dataFile<<"[";
    for (int memIdx=0;memIdx<numOfMemValues;memIdx++)
    {
        int bits = memorys[memIdx];
        logFile<<"memory: "<<bits<<" bits"<<endl;

        int registerNum=bits/6;

        for (int roundIdx = 0; roundIdx < round; roundIdx++){

            time_t seed = rand()+1316561;
            mt19937 rng(seed);
            uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
            hyperlogloglog::HyperLogLog<> hll (registerNum);
            for (key = 1; key <= nf; key++) {
                uint64_t y = dist(rng);
                hll.add(y);
            }

            double result= hll.estimate();
            values[memIdx][roundIdx]= result;

            if((roundIdx+1)%1==0){
                logFile<<"nf: "<<nf<<", round: "<<(roundIdx+1)<<", estimate: "<<result<<endl;
            }
        }
        double mean = stats_mean(values[memIdx], round);
        double coe=nf/mean;
        logFile<<"bits: "<<bits<<" nf: "<<nf<<", mean: "<<mean<<", coe: "<<coe<<endl;//实际基数
        coeSum+=coe;

        for (int roundIdx = 0; roundIdx < round; roundIdx++){
            values[memIdx][roundIdx] =values[memIdx][roundIdx]/nf;
        }

        double varer =stats_variance(values[memIdx],1.0, round);
        double stdErr=sqrt(varer);
        double MVP=varer*bits;
        logFile<<"bits: "<<bits<<" nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
        dataFile<<stdErr<<","<<flush;
    }
    dataFile<<"],"<<endl;

    dataFile<<"]"<<endl;
    logFile.close();
    dataFile.close();

    return 0;
}
