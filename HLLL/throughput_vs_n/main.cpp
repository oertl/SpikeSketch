#include <iostream>
#include <random>

#include "HyperLogLogLog.hpp"
#include <random>
#include <fstream>
#include "stats.h"
using namespace std;


//uint32_t calculateTotalMem(int registerNum,int tableSize){
//    double mem=3.0*registerNum+(ceil(std::log2(registerNum))+6)*tableSize+6.0;
//    return ceil(mem);
//}
//
//uint32_t calculateExpectedMem(int registerNum,double ratio){
//    double mem=(3.0+ratio)*registerNum+6.0;
//    return ceil(mem);
//}
//
////int findRegNum(uint32_t expectedMem,double ratio){
////    int lower=(expectedMem-6)/32;
////    int upper=(expectedMem-6)/2.5;
////
////    while(lower<upper){
////        int median=(lower+upper)/2;
////        uint32_t mem= calculateExpectedMem(median,ratio);
////        if(mem<expectedMem){
////            if(calculateExpectedMem(median+1,ratio)>expectedMem){
////                return median;
////            }
////            lower=median;
////        }else if(mem>expectedMem){
////            upper=median;
////        }else{
////            return median;
////        }
////    }
////    return lower;
////}
//
//int findRegNum(uint32_t expectedMem,uint32_t expectedTableSize){
//    int lower=(expectedMem-6)/32;
//    int upper=(expectedMem-6)/2.5;
//
//    while(lower<upper){
//        int median=(lower+upper)/2;
//        uint32_t mem= calculateTotalMem(median,expectedTableSize);
//        if(mem<expectedMem){
//            if(calculateTotalMem(median+1,expectedTableSize)>expectedMem){
//                return median;
//            }
//            lower=median;
//        }else if(mem>expectedMem){
//            upper=median;
//        }else{
//            return median;
//        }
//    }
//    return lower;
//}
//
////测试table size
//uint32_t test(int registerNum,uint32_t nf,int round){
//    uint64_t key;
//
//    double totalMaxSize=0.0;
//    for (int roundIdx = 0; roundIdx < round; roundIdx++){
//        time_t seed = rand()+1316561;
//        mt19937 rng(seed);
//        uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
//        hyperlogloglog::HyperLogLogLog<> hlll (registerNum);
//
//        double maxSize=0.0;
//        for (key = 1; key <= nf; key++) {
//            uint64_t y = dist(rng);
//            double size=hlll.add(y);
//            if(size>maxSize){
//                maxSize=size;
//            }
//        }
//        if(maxSize>totalMaxSize){
//            totalMaxSize=maxSize;
//        }
//    }
//    return totalMaxSize;
//}
////
//////获取RegNum
////int getRegNum(uint32_t memoryBits,uint32_t nf,int round){
////    cout<<"mem: "<<memoryBits<<", nf"<<nf<<", round: "<<round<<", get regNum"<<endl;
////    int regNumUpper= findRegNum(memoryBits,0.0);
////    uint32_t realTableSizeUpper=test(regNumUpper, nf,round);
////
////    double ratio=(ceil(std::log2(regNumUpper))+6.0)*realTableSizeUpper/regNumUpper;
////    ratio=ceil(ratio*100.0)/100.0;//保留两位小数
////
////    int regNum=findRegNum(memoryBits,ratio);
////    cout<<"regNum: "<<regNum<<endl;
////    return regNum;
////}
//
////获取RegNum
//int getRegNum(uint32_t memoryBits,uint32_t nf,int round){
//    cout<<"mem: "<<memoryBits<<", nf: "<<nf<<", round: "<<round<<", get regNum"<<endl;
//    int regNumUpper= findRegNum(memoryBits,0);
//    uint32_t realTableSizeUpper=test(regNumUpper, nf,round);
//
//    int regNumLower=findRegNum(memoryBits,realTableSizeUpper);
//    uint32_t realTableSizeLower=test(regNumLower, nf,round);
//
//    while(((double)realTableSizeUpper/realTableSizeLower>1.1)&&(realTableSizeUpper-realTableSizeLower>=1)){
//        regNumUpper=findRegNum(memoryBits,realTableSizeLower);
//        uint32_t newRealTableSizeUpper=test(regNumUpper, nf,round);
//        if(newRealTableSizeUpper>=realTableSizeUpper){
//            break;
//        }
//        if(newRealTableSizeUpper<=realTableSizeLower){
//            regNumLower=regNumUpper;
//            break;
//        }
//        realTableSizeUpper=newRealTableSizeUpper;
//
//        regNumLower=findRegNum(memoryBits,realTableSizeUpper);
//        uint32_t newRealTableSizeLower=test(regNumLower, nf,round);
//        if(newRealTableSizeLower<=realTableSizeLower){
//            break;
//        }
//        realTableSizeLower=newRealTableSizeLower;
//    }
//    cout<<"regNum: "<<regNumLower<<"table size: "<<realTableSizeUpper<<endl;
//    return regNumLower;
//}

int main() {
    srand(time(NULL));

    vector<uint64_t> nfs;


    for(int x=10;x<=25;x+=2){
        for(int i=0;i<4;i++){
            double temp=pow(2,x+0.5*i);
            if(temp>(1<<25)){
                break;
            }
            uint64_t nf=temp;
            nfs.push_back(nf);
        }
    }

    int numOftestCards=nfs.size();

    int memoryBits=1<<13;
    int roundTimes=100;

    uint64_t key;

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

//    int registerNum=getRegNum(memoryBits,100000,1000);
    int registerNum=2100;

    mt19937 rng(0);
    for (int nfi=0;nfi<numOftestCards;nfi++)
    {
        long double totalTimeNS=0.0;
        uint64_t nf = nfs[nfi];

        for(int roundIdx=0;roundIdx<roundTimes;roundIdx++) {
            __uint128_t resns;
            uint64_t* keys=new uint64_t [nf];
            uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
            for(int i=0;i<nf;i++){
                keys[i]=dist(rng);
            }
            hyperlogloglog::HyperLogLogLog<> hlll (registerNum);

            timespec time1, time2;
            clock_gettime(CLOCK_MONOTONIC, &time1);
            for (int k = 0; k < nf; k++) {
                hlll.add(keys[k]);
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
