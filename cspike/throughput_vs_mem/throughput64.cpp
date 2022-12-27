//
// Created by 李克剑 on 2022/7/12.
//
//
// Created by SamLee on 2022/6/23.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "spikesketch64.cpp"
#include <random>


int main(int argc,char* argv[]) {

    srand(time(NULL));

    int start=atoi(argv[1]);
    int end=atoi(argv[2]);

    int roundTimes=100;
    int nf=1<<20;
    uint64_t* keys=new uint64_t [nf];

    mt19937 rng(0);
    uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
    for(int i=0;i<nf;i++){
        keys[i]=dist(rng);
    }

    vector<int> memorys;
    for(int i=start;i<=end;i++){
        memorys.push_back(1<<i);
    }
    int numOfMemValues=memorys.size();
    double* values=new double [numOfMemValues];


    char traceFilePath[100];
    ofstream logFile;
    sprintf(traceFilePath, "%d-%d.log", start, end);
    logFile.open(traceFilePath,ios::out);

    ofstream dataFile;
    sprintf(traceFilePath, "%d-%d.dat", start, end);
    dataFile.open(traceFilePath,ios::out);

    logFile<<"nf: "<<nf<<endl;
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


    dataFile<<"[";
    for (int memIdx=0;memIdx<numOfMemValues;memIdx++)
    {
        int memoryBits = memorys[memIdx];

        int numOfBuckets=memoryBits>>6;
        int logNumOfBuckets=std::log2(numOfBuckets);

        logFile<<"memory: "<<memoryBits<<" bits, num of buckets: "<<numOfBuckets<<endl;
        long double totalTimeNS=0.0;

        for(int roundIdx=0;roundIdx<roundTimes;roundIdx++) {
            __uint128_t resns;

            time_t seed = rand()+1316561;
            spike_sketch ss = spike_sketch(logNumOfBuckets, seed);

            timespec time1, time2;
            clock_gettime(CLOCK_MONOTONIC, &time1);
            for (int k = 0; k < nf; k++) {
                ss.update(keys[k]);
            }
            clock_gettime(CLOCK_MONOTONIC, &time2);
            resns = (__uint128_t) (time2.tv_sec - time1.tv_sec) * 1000000000LL +
                    (time2.tv_nsec - time1.tv_nsec); //time spent on insert operation in nsec

            long double updateTimeNS=(long double)resns/nf;
            values[memIdx]=updateTimeNS;
            totalTimeNS+=updateTimeNS;
        }

        long double avgUpdateTimeNS=totalTimeNS/roundTimes;
        values[memIdx]=avgUpdateTimeNS;
        logFile<<"bits: "<<memoryBits<<" nf: "<<nf<<", avgUpdateTime(ns): "<<avgUpdateTimeNS<<endl;//实际基数
        dataFile<<avgUpdateTimeNS<<","<<flush;
    }
    dataFile<<"],"<<endl;

    dataFile<<"]"<<endl;
    logFile.close();
    dataFile.close();

    return 0;
}
