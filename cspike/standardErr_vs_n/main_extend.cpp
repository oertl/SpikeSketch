#include <iostream>
#include "utils/MurmurHash3.h"
#include "spike_sketch_extend.h"
#include <list>
#include <numeric>
#include "utils/stats.h"
#include <fstream>
#include <ctime>
#include <unordered_set>
#include <random>

using namespace std;

static mt19937 rng(time(NULL));//C++11的新随机数，可以生成int范围的随机数

double mutiBktQuery(vector<spike_sketch>& spikeSketchArray,double alpha0,double alpha1,double beta0,double beta1,double coe){
    double threshold=0.04;

    int numOfBuckets=spikeSketchArray.size();

    int numOfStages=0;
    int stageWidth=0;
    if(numOfBuckets>0){
        numOfStages=spikeSketchArray[0].getNumOfStages();
        stageWidth=spikeSketchArray[0].getStageWidth();
    }

    vector<double> emptyRatios;
    vector<int> numOfRegsInStages;
    int fzSum=0;

    for(int stageIdx=0;stageIdx<numOfStages;stageIdx++){
        double totalRegs=0;
        double unknownRegs[32] = {0};
        double rankRegs[32] = {0};
        double rankNonEmptyRegs[32] = {0};

        for(int bucketIdx=0;bucketIdx<numOfBuckets;bucketIdx++){
            spike_sketch ss=spikeSketchArray[bucketIdx];

            for(int cellIdx=stageIdx*stageWidth;cellIdx<(stageIdx+1)*stageWidth&&cellIdx<ss.n-1;cellIdx++){
                totalRegs++;
                if(!ss.tension(cellIdx)||ss.E[cellIdx]==1){
                    rankRegs[ss.S[cellIdx]] += 1;
                    if(ss.S[cellIdx]!=0){
                        rankNonEmptyRegs[ss.S[cellIdx]] += 1;
                    }
                }else{
                    if(ss.S[cellIdx]!=1){
                        unknownRegs[ss.S[cellIdx]] += 1;
                    }else{
                        rankRegs[ss.S[cellIdx]] += 1;
                    }
                }
            }
        }

        double smallerRegs = rankRegs[0] + rankRegs[1];
        double smallerNonEmpty = rankNonEmptyRegs[0] + rankNonEmptyRegs[1];
        double smallerNonEmptyRatio;
        for (int r = 2; r < 32; r++){
            smallerNonEmptyRatio = smallerNonEmpty / smallerRegs;
            smallerRegs += unknownRegs[r];
            smallerRegs += rankRegs[r];
            smallerNonEmpty += smallerNonEmptyRatio * unknownRegs[r];
            smallerNonEmpty += rankNonEmptyRegs[r];
        }

        double emptyRatio= 1- (double)smallerNonEmpty / smallerRegs;

        if(emptyRatio>=threshold){
            fzSum++;
            emptyRatios.push_back(emptyRatio);
        }
    }


    if(fzSum>2){//小流估计
        double sum=0;
        for(int stageIdx=0;stageIdx<numOfStages;stageIdx++){
            if(emptyRatios[stageIdx]>=threshold){
                double o_z = double(numOfStages-1-stageIdx)/numOfStages;
                sum+=-1.0*pow(4,o_z)*numOfBuckets*5*log(emptyRatios[stageIdx]);
            }
        }
        double result=((double)numOfStages/fzSum)*sum;
        return result;
    }else{
        double result=0;
        for(int bucketIdx=0;bucketIdx<numOfBuckets;bucketIdx++){
            spike_sketch ss=spikeSketchArray[bucketIdx];
            result+=ss.query(alpha0,alpha1,beta0,beta1,coe);
        }
        return result;
    }
}


int main() {

    srand(time(NULL));

    int round = 1000;

    vector<uint64_t> nfs;

    for(int x=10;x<15;x+=1){
        for(int i=0;i<2;i++){
            double temp=pow(2,x+0.5*i);
            uint64_t nf=temp;
            nfs.push_back(nf);
        }
    }

    int numOftestCards=nfs.size();

    int memoryBits=1<<13;
    int numOfBuckets=memoryBits>>6;

    uint64_t key;

    int n = 20; // need for a loop optimization
    int ncode = 4;
    int p = 12;


    double** values=new double*[numOftestCards];
    for(int i=0;i<numOftestCards;i++){
        values[i]=new double[round];
    }

    double alpha0=0.1;
    double alpha1=0.88;
    double beta0=1.12;
    double beta1=1.46;
    double myCoe=0.573;

    ofstream logFile;
    logFile.open("total.log",ios::out);

    ofstream dataFile;
    dataFile.open("total.dat",ios::out);


    logFile<<"ncodes: "<<ncode<<", ";
    logFile<<"p: "<<p<<", ";//概率p
    logFile<<"ncell: "<<n<<", ";//cell总个数，这里是二十
    logFile<<"round: "<<round<<endl;//实验次数
    logFile<<"alpha0: "<<alpha0<<", ";
    logFile<<"alpha1: "<<alpha1<<endl;
    logFile<<"beta0: "<<beta0<<", ";
    logFile<<"beta1: "<<beta1<<endl;
    logFile<<"beta1: "<<beta1<<endl;
    logFile<<"coe: "<<myCoe<<endl;
    logFile<<"mem: "<<memoryBits<<" bits"<<endl;
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
//        logFile<<"processing "<<nfi<<"cards: "<<nf<<endl;


        for (int roundIdx = 0; roundIdx < round; roundIdx++){
//            if((roundIdx+1)%(1000)==0){
//                logFile<<nf<<": "<<(roundIdx+1)<<"rounds"<<endl;
//            }

            time_t seed = rand()+1316561;

            vector<spike_sketch> spikeSketchArray;
            for(int bktIdx=0;bktIdx<numOfBuckets;bktIdx++){
                spike_sketch ss = spike_sketch(n, p, ncode, seed);
                spikeSketchArray.push_back(ss);
            }

            for (key = 1; key <= nf; key++) {
                uint32_t tempInt32=0;
                MurmurHash3_x86_32(&key, 4, seed+231321, &tempInt32);
                spike_sketch ss=spikeSketchArray[tempInt32%numOfBuckets];
                ss.update(key);

                if (!ss.valid()) {
                    for (int i = 0; i < ss.n; i++)
                        logFile << ss.S[i] << " ";
                    break;
                }
            }

            double result= mutiBktQuery(spikeSketchArray,alpha0,alpha1,beta0,beta1,myCoe);
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
        double MVP=varer*64.0*numOfBuckets;
        logFile<<"nf: "<<nf<<", stdErr: "<<stdErr<<", MVP: "<<MVP<<endl;//实际基数
        dataFile<<stdErr<<","<<flush;
    }
    dataFile<<"],"<<endl;
    dataFile<<"]"<<endl;

    logFile.close();
    dataFile.close();

    return 0;
}
