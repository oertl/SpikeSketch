#include <iostream>

#include "cpc_sketch.hpp"
#include "stats.h"
#include <fstream>

using namespace std;

uint32_t calculateTotalMem(int registerNum,int tableSize){
    double mem=8.0*registerNum+(ceil(std::log2(registerNum))+6)*tableSize+6.0;
    return ceil(mem);
}

int findRegNum(uint32_t expectedMem,uint32_t expectedTableSize){
    int lower=(expectedMem-6)/32;
    int upper=(expectedMem-6)/8;

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

        datasketches::cpc_sketch_alloc<std::allocator<uint8_t>> cpc_sk(registerNum);

        double maxSize=0.0;
        for (key = 1; key <= nf; key++) {

            uint64_t y = dist(rng);
            uint8_t lg_size=cpc_sk.update(y);
            double size=(double)(1<<lg_size);
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


//获取RegNum
int getRegNum(uint32_t memoryBits,uint32_t nf,int round){
    cout<<"mem: "<<memoryBits<<", nf"<<nf<<", round: "<<round<<", get regNum"<<endl;
    uint32_t iniTableSize=4;
    int regNumUpper= findRegNum(memoryBits,iniTableSize);//预设的表越小，reg越多
    uint32_t realTableSizeUpper=test(regNumUpper, nf,round);//reg越多，table size越大（猜测是这样）

//    cout<<realTableSizeUpper<<endl;
    cout<<"regNum1: "<<regNumUpper<<", table size1: "<<realTableSizeUpper<<endl;

    int regNumLower= findRegNum(memoryBits,realTableSizeUpper);
    uint32_t realTableSizeLower=test(regNumLower, nf,round);

    cout<<"regNum2: "<<regNumLower<<", table size2: "<<realTableSizeLower<<endl;

    if(realTableSizeLower>=realTableSizeUpper){
        return regNumLower;
    }else{

        regNumUpper=findRegNum(memoryBits,realTableSizeLower);
        realTableSizeUpper=test(regNumUpper, nf,round);
        cout<<"regNum3: "<<regNumUpper<<", table size3: "<<realTableSizeUpper<<endl;

        if(realTableSizeLower>=realTableSizeUpper) {
            return regNumUpper;
        }else{
            cout<<"memoryBits: "<<memoryBits<<", nf: "<<nf<<", round: "<<round<<" trouble encountered!"<<endl;

            regNumLower= findRegNum(memoryBits,realTableSizeUpper);
            regNumUpper= findRegNum(memoryBits,realTableSizeLower);
            while(regNumLower<regNumUpper){
                int testRegNum=(regNumLower+regNumUpper)/2;

                uint32_t realTableSize=test(testRegNum, nf,round);
                if(realTableSize>realTableSizeLower){
                    regNumUpper=testRegNum;
                }else{
                    cout<<"regNum4: "<<testRegNum<<", table size4: "<<realTableSize<<endl;
                    return testRegNum;
                }
            }
            return regNumLower;
        }
    }
}


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


    int registerNums[]={3,6,15,39,75,199,391,767,1503,2904,5759};

    dataFile<<"["<<endl;
    for (int memIdx=0;memIdx<numOfMemValues;memIdx++)
    {
        int memoryBits = memorys[memIdx];

//        int registerNum=getRegNum(memoryBits,100000,10000);
        int registerNum=registerNums[memIdx];

        logFile<<"memory: "<<memoryBits<<" bits, register num: "<<registerNum<<endl;

        long double totalTimeNS=0.0;
        for(int roundIdx=0;roundIdx<roundTimes;roundIdx++) {
            __uint128_t resns;

            datasketches::cpc_sketch_alloc<std::allocator<uint8_t>> cpc_sk(registerNum);

            timespec time1, time2;
            clock_gettime(CLOCK_MONOTONIC, &time1);
            for (int k = 0; k < nf; k++) {
                cpc_sk.update(keys[k]);
            }
            clock_gettime(CLOCK_MONOTONIC, &time2);
            resns = (__uint128_t) (time2.tv_sec - time1.tv_sec) * 1000000000LL +
                    (time2.tv_nsec - time1.tv_nsec); //time spent on insert operation in nsec

            long double updateTimeNS=(long double)resns/nf;
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
