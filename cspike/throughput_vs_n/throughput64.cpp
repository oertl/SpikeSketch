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


int main() {

    srand(time(NULL));

    vector<uint64_t> nfs;

    for(int x=10;x<25;x+=1){
        for(int i=0;i<2;i++){
            double temp=pow(2,x+0.5*i);
            uint64_t nf=temp;
            nfs.push_back(nf);
        }
    }
    nfs.push_back((1<<25));
    int numOftestCards=nfs.size();

    int memoryBits=1<<13;
    int roundTimes=100;

    double* values=new double [numOftestCards];

    ofstream logFile;
    logFile.open("total.log",ios::out);

    ofstream dataFile;
    dataFile.open("total.dat",ios::out);


    logFile<<"expected total Mem: "<<memoryBits<<"bits"<<endl;
    logFile<<"######################################\n\n\n"<<endl;


    dataFile<<"[";
    dataFile<<endl;

    dataFile<<"[";
    for(int nfi=0;nfi<numOftestCards;nfi++){
        dataFile<<nfs[nfi]<<",";
    }
    dataFile<<"],"<<endl;

    dataFile<<"["<<endl;

    int numOfBuckets=memoryBits>>6;
    int logNumOfBuckets=std::log2(numOfBuckets);

    mt19937 rng(0);
    for (int nfi=0;nfi<numOftestCards;nfi++)
    {
        uint64_t nf = nfs[nfi];
        long double totalTimeNS=0.0;


        for(int roundIdx=0;roundIdx<roundTimes;roundIdx++) {
            __uint128_t resns;

            uint64_t* keys=new uint64_t [nf];
            uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
            for(int i=0;i<nf;i++){
                keys[i]=dist(rng);
            }
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
            totalTimeNS+=updateTimeNS;
        }
        long double avgUpdateTimeNS=totalTimeNS/roundTimes;
        values[nfi]=avgUpdateTimeNS;

        logFile<<"bits: "<<memoryBits<<" nf: "<<nf<<", avgUpdateTime(ns): "<<avgUpdateTimeNS<<endl;//实际基数
        dataFile<<avgUpdateTimeNS<<","<<flush;
    }
    dataFile<<"],"<<endl;

    dataFile<<"]"<<endl;
    logFile.close();
    dataFile.close();

    return 0;
}

