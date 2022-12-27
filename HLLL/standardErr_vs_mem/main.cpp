#include <iostream>
#include <random>

#include "HyperLogLogLog.hpp"
#include <random>
#include <fstream>
#include "stats.h"
using namespace std;


uint32_t calculateTotalMem(int registerNum,int tableSize){
    double mem=3.0*registerNum+(ceil(std::log2(registerNum))+6)*tableSize+6.0;
    return ceil(mem);
}

uint32_t calculateExpectedMem(int registerNum,double ratio){
    double mem=(3.0+ratio)*registerNum+6.0;
    return ceil(mem);
}

//int findRegNum(uint32_t expectedMem,double ratio){
//    int lower=(expectedMem-6)/32;
//    int upper=(expectedMem-6)/2.5;
//
//    while(lower<upper){
//        int median=(lower+upper)/2;
//        uint32_t mem= calculateExpectedMem(median,ratio);
//        if(mem<expectedMem){
//            if(calculateExpectedMem(median+1,ratio)>expectedMem){
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

int findRegNum(uint32_t expectedMem,uint32_t expectedTableSize){
    int lower=(expectedMem-6)/32;
    int upper=(expectedMem-6)/2.5;

    while(lower<upper){
        int median=(lower+upper)/2;
        uint32_t mem= calculateTotalMem(median,expectedTableSize);
        if(mem<expectedMem){
            if(calculateTotalMem(median+1,expectedTableSize)>expectedMem){
                return median;
            }
            lower=median;
        }else if(mem>expectedMem){
            upper=median;
        }else{
            return median;
        }
    }
    return lower;
}

//测试table size
uint32_t test(int registerNum,uint32_t nf,int round){
    uint64_t key;

    double totalMaxSize=0.0;
    for (int roundIdx = 0; roundIdx < round; roundIdx++){
        time_t seed = rand()+1316561;
        mt19937 rng(seed);
        uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
        hyperlogloglog::HyperLogLogLog<> hlll (registerNum);

        double maxSize=0.0;
        for (key = 1; key <= nf; key++) {
            uint64_t y = dist(rng);
            double size=hlll.add(y);
            if(size>maxSize){
                maxSize=size;
            }
        }
        if(maxSize>totalMaxSize){
            totalMaxSize=maxSize;
        }
    }
    return totalMaxSize;
}
//
////获取RegNum
//int getRegNum(uint32_t memoryBits,uint32_t nf,int round){
//    cout<<"mem: "<<memoryBits<<", nf"<<nf<<", round: "<<round<<", get regNum"<<endl;
//    int regNumUpper= findRegNum(memoryBits,0.0);
//    uint32_t realTableSizeUpper=test(regNumUpper, nf,round);
//
//    double ratio=(ceil(std::log2(regNumUpper))+6.0)*realTableSizeUpper/regNumUpper;
//    ratio=ceil(ratio*100.0)/100.0;//保留两位小数
//
//    int regNum=findRegNum(memoryBits,ratio);
//    cout<<"regNum: "<<regNum<<endl;
//    return regNum;
//}

//获取RegNum
int getRegNum(uint32_t memoryBits,uint32_t nf,int round){
    cout<<"mem: "<<memoryBits<<", nf: "<<nf<<", round: "<<round<<", get regNum"<<endl;
    int regNumUpper= findRegNum(memoryBits,0);
    uint32_t realTableSizeUpper=test(regNumUpper, nf,round);

    int regNumLower=findRegNum(memoryBits,realTableSizeUpper);
    uint32_t realTableSizeLower=test(regNumLower, nf,round);

    while(((double)realTableSizeUpper/realTableSizeLower>1.1)&&(realTableSizeUpper-realTableSizeLower>=1)){
        regNumUpper=findRegNum(memoryBits,realTableSizeLower);
        uint32_t newRealTableSizeUpper=test(regNumUpper, nf,round);
        if(newRealTableSizeUpper>=realTableSizeUpper){
            break;
        }
        if(newRealTableSizeUpper<=realTableSizeLower){
            regNumLower=regNumUpper;
            break;
        }
        realTableSizeUpper=newRealTableSizeUpper;

        regNumLower=findRegNum(memoryBits,realTableSizeUpper);
        uint32_t newRealTableSizeLower=test(regNumLower, nf,round);
        if(newRealTableSizeLower<=realTableSizeLower){
            break;
        }
        realTableSizeLower=newRealTableSizeLower;
    }
    cout<<"regNum: "<<regNumLower<<"table size: "<<realTableSizeUpper<<endl;
    return regNumLower;
}

int main(int argc,char* argv[]) {
    int start=atoi(argv[1]);
    int end=atoi(argv[2]);


    srand(time(NULL));

    int round = 1000;

    vector<int> memorys;

    for(int i=start;i<=end;i++){
        memorys.push_back(1<<i);
    }

    uint64_t key;
    int nf=1<<20;


    int numOfMemValues=memorys.size();
    double** values=new double*[numOfMemValues];
    for(int i=0;i<numOfMemValues;i++){
        values[i]=new double[round];
    }


    char traceFilePath[100];
    ofstream logFile;
    sprintf(traceFilePath, "%d-%d.log", start, end);
    logFile.open(traceFilePath,ios::out);

    ofstream dataFile;
    sprintf(traceFilePath, "%d-%d.dat", start, end);
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

    dataFile<<"["<<endl;


    for (int memIdx=0;memIdx<numOfMemValues;memIdx++)
    {
        int memoryBits = memorys[memIdx];

        double avgSize=0.0;
        double totalMaxSize=0.0;
        int registerNum=getRegNum(memoryBits,100000,10000);
        logFile<<"memory: "<<memoryBits<<" bits, register num: "<<registerNum<<endl;
        for (int roundIdx = 0; roundIdx < round; roundIdx++){

            time_t seed = rand()+1316561;
            mt19937 rng(seed);
            uniform_int_distribution<uint64_t> dist(0,UINT64_MAX);
            hyperlogloglog::HyperLogLogLog<> hlll (registerNum);

            double maxSize=0.0;
            for (key = 1; key <= nf; key++) {

                uint64_t y = dist(rng);
                double size=hlll.add(y);
                if(size>maxSize){
                    maxSize=size;
                }
            }

            avgSize+= maxSize;
            if(maxSize>totalMaxSize){
                totalMaxSize=maxSize;
            }

            double result= hlll.estimate();
            values[memIdx][roundIdx]= result;
        }
        avgSize=avgSize/ round;

        double mean = stats_mean(values[memIdx], round);
        double coe=nf/mean;

        double totalAvgMem=3.0*registerNum+(ceil(std::log2(registerNum))+6)*avgSize+6;
        double totalMaxMem=3.0*registerNum+(ceil(std::log2(registerNum))+6)*totalMaxSize+6;

        logFile<<"nf: "<<", mean: "<<mean<<", coe: "<<coe<<", regsiterNum: "<<registerNum<<", avgSize: "<<avgSize<<", maxSize: "<<totalMaxSize<<", totalAvgMem: "<<totalAvgMem<<", totalMaxMem: "<<totalMaxMem<<endl;//实际基数
        coeSum+=coe;

        for (int roundIdx = 0; roundIdx < round; roundIdx++){
            values[memIdx][roundIdx] =values[memIdx][roundIdx]/nf;
        }

        double varer =stats_variance(values[memIdx],1.0, round);
        double stdErr=sqrt(varer);
        double MVP=varer*memoryBits;
        logFile<<"bits: "<<memoryBits<<" nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
        dataFile<<stdErr<<","<<flush;
    }
    dataFile<<"],"<<endl;
    dataFile<<"]"<<endl;

    logFile.close();
    dataFile.close();

    return 0;
}
