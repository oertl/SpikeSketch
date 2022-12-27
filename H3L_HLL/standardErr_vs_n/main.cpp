#include <iostream>
#include <random>

#include "HyperLogLog.hpp"
#include <random>
#include <fstream>
#include "stats.h"
using namespace std;

int main() {
    srand(time(NULL));

    int round = 1000;

    vector<uint64_t> nfs;

    for(int x=10;x<13;x+=1){
        for(int i=0;i<2;i++){
            double temp=pow(2,x+0.5*i);
            uint64_t nf=temp;
            nfs.push_back(nf);
        }
    }

    int numOftestCards=nfs.size();//sizeof(nfs)/sizeof(nfs[0])


    int memoryBits=1<<13;
    int registerNum=memoryBits/6;

    uint64_t key;

    double** values=new double*[numOftestCards];
    for(int i=0;i<numOftestCards;i++){
        values[i]=new double[round];
    }

    ofstream logFile;
    logFile.open("total.log",ios::out);

    ofstream dataFile;
    dataFile.open("total.dat",ios::out);


    logFile<<"registerNum: "<<registerNum<<", ";
    logFile<<"round: "<<round<<endl;//实验次数
    logFile<<"######################################\n\n\n"<<endl;


    dataFile<<"[";
    dataFile<<endl;

    dataFile<<"[";
    for(int nfi=0;nfi<numOftestCards;nfi++){
        dataFile<<nfs[nfi]<<",";
    }
    dataFile<<"],"<<endl;

    double coeSum=0.0;

    dataFile<<"["<<endl;
    for (int nfi=0;nfi<numOftestCards;nfi++)
    {
        uint64_t nf = nfs[nfi];

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
            values[nfi][roundIdx]= result;
        }

        double mean = stats_mean(values[nfi], round);
        double coe=nf/mean;
        logFile<<"nf: "<<nf<<", mean: "<<mean<<", coe: "<<coe<<endl;//实际基数
        coeSum+=coe;

        for (int roundIdx = 0; roundIdx < round; roundIdx++){
            values[nfi][roundIdx] =values[nfi][roundIdx]/nf;
        }
        double varer =stats_variance(values[nfi],1.0, round);
        double stdErr=sqrt(varer);
        double MVP=varer*6.0*registerNum;
        logFile<<"nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
        dataFile<<stdErr<<","<<flush;
    }
    dataFile<<"],"<<endl;
    dataFile<<"]"<<endl;

    logFile.close();
    dataFile.close();

    return 0;
}
