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


int main() {
    srand(time(NULL));

    int round = 1000;

    uint64_t nf=1<<25;

    vector<int> regNums={49,193,1201,4803};
    int numOfRegNums=regNums.size();


    uint64_t key;

    double** values=new double*[numOfRegNums];
    for(int i=0;i<numOfRegNums;i++){
        values[i]=new double[round];
    }

    ofstream logFile;
    logFile.open("total.log",ios::out);

    ofstream dataFile;
    dataFile.open("total.dat",ios::out);


    logFile<<"nf: "<<nf<<endl;
    logFile<<"round: "<<round<<endl;//实验次数
    logFile<<"######################################\n\n\n"<<endl;


    dataFile<<"[";
    dataFile<<endl;

    dataFile<<"[";
    for(int i=0;i<numOfRegNums;i++){
        dataFile<<regNums[i]<<",";
    }
    dataFile<<"],"<<endl;


    double coeSum=0.0;

    dataFile<<"["<<endl;
    for (int i=0;i<numOfRegNums;i++)
    {
        double avgSize=0.0;
        double totalMaxSize=0.0;
        int registerNum=regNums[i];
        for (int roundIdx = 0; roundIdx < round; roundIdx++){

            if((roundIdx+1)%100==0){
                logFile<<"nf: "<<nf<<", round: "<<(roundIdx+1)<<endl;
            }

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

            avgSize+= maxSize;
            if(maxSize>totalMaxSize){
                totalMaxSize=maxSize;
            }

            double result= cpc_sk.get_estimate();
            values[i][roundIdx]= result;
        }
        avgSize=avgSize/ round;


        double mean = stats_mean(values[i], round);
        double coe=nf/mean;

        double totalAvgMem=8*registerNum+(ceil(std::log2(registerNum))+6)*avgSize+6;
        double totalMaxMem=8*registerNum+(ceil(std::log2(registerNum))+6)*totalMaxSize+6;

        logFile<<"nf: "<<", mean: "<<mean<<", coe: "<<coe<<", regsiterNum: "<<registerNum<<", avgSize: "<<avgSize<<", maxSize: "<<totalMaxSize<<", totalAvgMem: "<<totalAvgMem<<", totalMaxMem: "<<totalMaxMem<<endl;//实际基数
        coeSum+=coe;


        for (int roundIdx = 0; roundIdx < round; roundIdx++){
            values[i][roundIdx] =values[i][roundIdx]/nf;
        }
        double varer =stats_variance(values[i],1.0, round);
        double stdErr=sqrt(varer);
        double MVP=varer*totalMaxMem;
        logFile<<"nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
        dataFile<<stdErr<<","<<flush;
    }
    dataFile<<"],"<<endl;
    dataFile<<"]"<<endl;

    logFile.close();
    dataFile.close();

    return 0;
}
