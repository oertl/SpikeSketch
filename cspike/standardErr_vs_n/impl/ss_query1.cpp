//
// Created by SamLee on 2022/6/15.
//
#include "../spike_sketch_extend.h"
#include <cmath>       /* log2 */
#include "../utils/stats.h"
int nstage;
double stage_width;
double *stage_p;

spike_sketch::spike_sketch(int n, int p, int ncodes, uint32_t seed){
    this->n = n;//cell数量，default=41
    this->p = p;//更新概率，default=12
    this->q = double(16)/double(16-this->p);
    this->ncodes = ncodes;//如果全为一级编码，可以表示的取值数量.log2nocdes为单个cell的位长，default=4
    this->seed = seed;//随机数种子
    this->seed2 = seed^0xf0f0f0f0f;//随机数种子2

    nstage = log2(this->q)*2;
    stage_width = double(n)/nstage;
    stage_p = new double [nstage];
    for (int i = 0; i < nstage; i ++){
        stage_p[i] = 0xffffffff * pow(this->q, -double(nstage-i-1)/nstage);
    }


    this->S = new int [n];
    this->E = new int [n];
//    this->T = new bool [n];
    for (int i = 0; i < n; i ++){
        this->S[i] = 0;
        this->E[i] = 1;//如果张紧了，表示是否取到了当前值的位置；如果未张紧，表示是否取到当前值减一的位置。
//        this->T[i] = false;//对应索引的cell是否张紧
    }
    this->vlimit = int((ncodes-2) / 2);
    this->slimit = this->vlimit + ncodes;
}


bool spike_sketch::valid() {
//    for(int i=0;i<this->n;i++)
//        std::cout<<this->S[i]<<" ";
    int i = 0;
    double delta_1;
    for(int i=0;i<this->n-1;i++){
        delta_1=this->S[(i+1)%this->n]-this->S[i];
        if(abs(delta_1)>this->slimit)
            return false;
        else if (delta_1>this->vlimit){
            if (delta_1+this->S[(i+2)%this->n]-this->S[(i+1)%this->n]!=0)
                return false;
        }
        else if(-delta_1>this->vlimit){
            int left_i=i==0?n-1:i-1;
            if (delta_1+this->S[i]-this->S[left_i]!=0)
                return false;
        }
    }
    return true;
}

int spike_sketch::rho(uint64_t x) {
    if(x==0) return 31;
    uint32_t temp=(__builtin_clzll(x))/2;
    if(temp>=31) return 31;

    return  temp+ 1;
}

bool spike_sketch::tension(int j) {
    int r_delta, l_delta;
    r_delta = j != this->n-1 ? this->S[j] - this->S[j+1] : this->S[j] - this->S[0];
    l_delta = j != 0 ? this->S[j] - this->S[j-1]: this->S[j] - this->S[this->n-1];
    if (-r_delta >= this->vlimit or -l_delta >= this->vlimit)
        return true;
    else
        return false;
}

void spike_sketch::update(int key) {
    uint32_t x2;
    MurmurHash3_x86_32(&key, 4, this->seed2, &x2);
    int j;//cell索引
    double v=0;//更新值
    j = x2 % this->n;

    int stage_i = int(j / stage_width);

    if (x2 < stage_p[stage_i]){
        uint64_t x[2];
        MurmurHash3_x64_128(&key, 8, this->seed, &x);
        v += this->rho(x[0]);
    }

    if (v == this->S[j] and this->E[j] == 0 and tension(j))
        this->E[j] = 1;
    if (v == this->S[j]-1 and this->E[j] == 0 and !tension(j))
        this->E[j] = 1;

    if (v > this->S[j]){
//        cout<<"j: "<<j<<" v: "<<v<<" key: "<<key<<"\n";
        this->adjust(j, v);
//        this->updateT();
//        for(int i=0;i<this->n;i++)
//            std::cout<<this->S[i]<<" ";
//        std::cout<<endl;
//        for(int i=0;i<this->n;i++)
//            std::cout<<this->E[i]<<" ";
//        std::cout<<endl;
//        for(int i=0;i<this->n;i++)
//            std::cout << (!this->T[i] ? "F" : "T") << " ";
//        std::cout<<endl;
    }
}

//void spike_sketch::updateT() {
//    for (int i = 0; i < this->n ; i++)
//        this->T[i] = tension(i);
//}

void spike_sketch::extended_inc(int j, int v)
{
    int old_j=this->S[j];
    bool old_tension=this->tension(j);
    this->S[j] = v;
    if (v<old_j){
        cout<<"what?!";
        exit(1);
    }
    if (v > old_j){
        if (tension(j)){
            this->E[j] = 0;
        }
        else{
            if(v-old_j==1 and  (!old_tension or this->E[j] == 1)){
                this->E[j] = 1;
            }
            else{
                this->E[j] = 0;
            }
        }
    }
}

// each word is 64 bits, q=4, ncode=4, tension = 1
// use a loop transform, the last delta can be inferred or using the other spike, thus, it only requires 1bit to encode
// the query function with tension can be improved
void spike_sketch::adjust(int j, int v) {
    int old_j=this->S[j];
    bool old_tension=this->tension(j);
    this->S[j] = v;
    if (tension(j)){
        this->E[j] = 1;
    }
    else{
        if(v-old_j==1 and (!old_tension or this->E[j] == 1)){
            this->E[j] = 1;
        }
        else{
            this->E[j] = 0;
        }
    }

    int k;
    int delta_1, delta_2, min_delta;
    int start_k = 1;
    if (j==0){
        min_delta = std::min(v - this->S[this->n-1], v - this->S[j+1]);
        if (min_delta > this->vlimit){
            this->extended_inc(this->n-1, v-min(min_delta, this->slimit));
            this->extended_inc(j+1, v-min(min_delta, this->slimit));
            start_k = 2;
        }
    }
    else if(j==this->n-1){
        min_delta = std::min(v - this->S[j-1], v - this->S[0]);
        if (min_delta > this->vlimit){
            this->extended_inc(j-1, v-min(min_delta, this->slimit));
            this->extended_inc(0, v-min(min_delta, this->slimit));
            start_k = 2;
        }
    }
    else{
        min_delta = std::min(v - this->S[j-1], v - this->S[j+1]);
        if (min_delta > this->vlimit){
            this->extended_inc(j-1, v-min(min_delta, this->slimit));
            this->extended_inc(j+1, v-min(min_delta, this->slimit));
            start_k = 2;
        }
    }
    k = start_k;
    while (k < this->n-1){
        delta_1 = this->S[(j+k)%this->n] - this->S[(j+k-1)%this->n];
        delta_2 = this->S[(j+k+1)%this->n] - this->S[(j+k)%this->n];

        if (delta_1 > this->vlimit){
            if (delta_1 + delta_2 > 0){
                std::cout<<"error!"<<std::endl;
                cout<<"j: "<<j<<" v: "<<v<<" seed: "<<this->seed<<"\n";
                for(int i=0;i<this->n;i++)
                    std::cout<<this->S[i]<<" ";
                std::cout<<endl;
                exit(1);
            }
            else if (delta_1 + delta_2 < 0){
                this->extended_inc((j+k+1)%this->n, this->S[(j+k-1)%this->n]);
                k+=2;
                continue;
            }
            else{
                break;
            }
        }
        else{
            if (delta_1 < -this->vlimit){
                this->extended_inc((j+k)%this->n, this->S[(j+k-1)%this->n] - this->vlimit);
                k += 1;
                continue;
            }
            else{
                if (abs(delta_2) <= this->vlimit)
                    break;
                else{
                    k+=1;
                    continue;
                }
            }
        }
    }

    k = start_k;
    while (k < this->n-1){
        int current_i=j-k<0?j-k+n:j-k;
        int left_i=j-k-1<0?j-k-1+n:j-k-1;
        int right_i=j-k+1<0?j-k+1+n:j-k+1;
        delta_1 = this->S[current_i] - this->S[right_i];
        delta_2 = this->S[left_i] - this->S[current_i];

        if (delta_1 > this->vlimit){
            if (delta_1 + delta_2 > 0){
                std::cout<<"error!"<<std::endl;
                cout<<"j: "<<j<<" v: "<<v<<" seed: "<<this->seed<<"\n";
                for(int i=0;i<this->n;i++)
                    std::cout<<this->S[i]<<" ";
                std::cout<<endl;
                exit(1);
            }
            else if(delta_1 + delta_2 < 0){
                this->extended_inc(left_i, this->S[right_i]);
                k+=2;
                continue;
            }
            else{
                break;
            }
        }
        else{
            if (delta_1 < -this->vlimit){
                this->extended_inc(current_i, this->S[right_i] - this->vlimit);
                k += 1;
                continue;
            }
            else{
                if (abs(delta_2) <= this->vlimit)
                    break;
                else{
                    k+=1;
                    continue;
                }
            }
        }
    }
}

double spike_sketch::query(double alpha0,double alpha1,double beta0,double beta1,double coe) {
    double sum=0;

    for (int i = 0; i < this->n; i++) {
        sum += 1.0/getF_R_j(i,alpha0,alpha1,beta0,beta1);
    }

    return coe*(this->n*this->n/sum);
}

double spike_sketch::getF_R_j(int j,double alpha0,double alpha1,double beta0,double beta1){
    int d_=this->n-1;//d'

    int R_j=this->S[j];

    int stagej = int(j / stage_width);//the index of stage
    double o_j = double(nstage-1-stagej)/nstage;
    int bit=this->E[j];

    if(!tension(j) && j<d_){
        return pow(this->q,o_j)*(bit* getZ_k_1(R_j,alpha1)+(1-bit)* getZ_k_0(R_j,alpha0));
    }
    if(tension(j) && j<d_ && (bit==1)){
        return pow(this->q,o_j)*getZ_k(R_j);
    }
    if(tension(j) && j<d_ && (bit==0)){
        return pow(this->q,o_j)*getZ_k_minus(R_j,alpha0,alpha1);
    }
    
    if(!tension(j) && j>=d_){
        return pow(this->q,o_j)*beta1*getZ_k(R_j);
    }
    if(tension(j) && j>=d_){
        return pow(this->q,o_j)*beta0*getZ_k(R_j);
    }
}

double spike_sketch::getZ_k_0(double k,double alpha0) {
    return alpha0*pow(4.0,k);
}

double spike_sketch::getZ_k_1(double k,double alpha1){
    return alpha1*pow(4.0,k);
}

double spike_sketch::getZ_k(double k) {
    return pow(4.0,k);
}


// 只考虑 前d'的rank平均
double spike_sketch::getZ_k_minus(double k,double alpha0,double alpha1) {
    int d_=this->n-1;//d'

    double sum=0.0;
    int num=0;

    for (int i = 0; i < d_; i++){
        if(this->S[i]<k && (!tension(i)||this->E[i]==1)){
            sum+=this->S[i];
            num++;
        }
    }

    if (num==0){
        return getZ_k(k-1);
    }else{
        double mean=sum/(double)num;
        return pow(4.0,mean);
    }
}



int spike_sketch::getNumOfStages(){
    return nstage;
}

int spike_sketch::getStageWidth(){
    return stage_width;
}