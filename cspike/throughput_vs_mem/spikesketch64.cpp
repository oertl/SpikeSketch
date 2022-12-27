#include <iostream>
#include <cmath>       /* log2 */
#include "MurmurHash2.h"
using namespace std;
#define NSTAGE 4
#define N 20
#define STAGE_WIDTH 5
#define VLIMIT 1
#define SLIMIT 5
#define Q 4
/*
00:-1
01:0
10:1
1100:2,-2
1101:3,-3
1110:4,-4
1111:5,-5
*/
class spike_sketch{
    int *stage_p;//各stage对应的采样率
    int m;
    int e;
    uint32_t seed;
    uint32_t seed2;
    uint64_t *S;
    unsigned char* Sp;//sketch指针
    unsigned char* Se;
    unsigned char* Sbase;
public:
    spike_sketch(int e,uint32_t seed){
        this->seed=seed;
        this->seed2 = seed^0xf0f0f0f0f;
        this->e=e;
        this->m= pow(2,e);
        /*5bit base + 1*41bit extension + 2*41bit cell */
//        this->S=0;
        /*将所有的E置为1*/
//        *(this->Sp+10)=0b11111100;
//        *(uint32_t*)(this->Sp+11)=0xFFFFFFFF;
//        *(this->Sp+15)=0b00000111;
        S=new uint64_t[m];
        for(int i=0;i<m;i++){
            this->S[i]=0x07FFFF5555555555;
        }
        this->stage_p = new int [NSTAGE];
        stage_p[0]=256;
        stage_p[1]=181;
        stage_p[2]=128;
        stage_p[3]=90;
    }
    void update(uint64_t key){
        uint32_t x2=MurmurHash2(&key, 8, this->seed2);
        unsigned int i,o,j;
        i=x2&(this->m-1);
        o=(x2>>this->e)&0xFF;
        j = (x2 >>(e+8))%N;

        int stage_i = int(j / STAGE_WIDTH);
        if (o >= this->stage_p[stage_i])
            return;


//        unsigned __int128 x;
//        MurmurHash3_x64_128(key, 4, this->seed, &x);
        Sp=(unsigned char*)(S+i);
        Se=Sp+5;
        Sbase=Sp+7;
        uint64_t x=MurmurHash64A(&key, 8, this->seed);
        double v= this->rho(x);


        int Sj=get_cell(j);
//        bool Ej=get_extension(j);
        bool Tj=tension(j);
        if (v==(Sj-(!Tj))){
            set_extension(j,1);
        }
//        if (v == Sj and !Ej and Tj)
//            set_extension(j,1);
//        if (v == Sj-1 and !Ej and !Tj)
//            set_extension(j,1);
//        cout<<"j: "<<j<<" v: "<<v<<" key: "<<key<<"\n";
        if (v > Sj){
//            cout<<"j: "<<j<<" v: "<<v<<"\n";
            this->adjust(j, v-Sj);
//            for(int i=0;i<N;i++)
//                std::cout<<this->get_cell(i)<<" ";
//            std::cout<<endl;
//        for(int i=0;i<N;i++)
//            std::cout<<this->get_extension(i)<<" ";
//        std::cout<<endl;
        }
    }
    double query(){
        double x=0;
        int stagei;
        x += 1/pow(Q, get_cell(N-1));//最后两个cell无tesnsion位
        x += 1/pow(Q, get_cell(N-2));
        for (int i = 0; i < N-2; i++){
            stagei = int(i / STAGE_WIDTH);
            int Si=get_cell(i);
            if (!tension(i)){
                x += 1/pow(Q,  Si+ double(NSTAGE-1-stagei)/NSTAGE);
                x += (1-(this->get_extension(i)?1:0)) * 1/pow(Q, Si + double(NSTAGE-1-stagei)/NSTAGE - 1);
            }
            else
            {
                if (get_extension(i))
                    x += 1/pow(Q, Si + double(NSTAGE-1-stagei)/NSTAGE);
                else{
                    double mean_x=0;
                    int xs=0;
                    for(int l =0;l<N;l++) {
                        int Sl=get_cell(l);
                        if(Sl < Si and (!tension(l) or get_extension(l))){
                            mean_x+=Sl;
                            xs++;
                        }
                    }
                    x += 1/pow(Q, xs==0?(Si + double(NSTAGE-1-stagei)/NSTAGE - 1):(mean_x/xs+ double(NSTAGE-1-stagei)/NSTAGE));
                }

            }
        }
        return (N*N/x)*0.558*0.995;
    }
    int get_cell(int j){//0<=j<=40
//        int v=(int) (this->S>>123);//base
        int v=(int) (*(Sbase))>>3;
        if (j==0){
            return v;
        }
        int i;
        if (j<N/2) {
            i=0;
            if (is_spike_code(i)) {//第一个delta可能是spikecode的后半部分。
                v-= (get_code(0)+2);
//                switch (get_code(0)) {
//                    case '\000':
//                        v -= 2;
//                        break;
//                    case '\001':
//                        v -= 3;
//                        break;
//                    case '\002':
//                        v -= 4;
//                        break;
//                    case '\003':
//                        v -= 5;
//                        break;
//                }
                i++;
            }
            while (i < j) {
                char codei= get_code(i);
                if(codei=='\003'){
                    if(i==j-1){
                        v+=(get_code(j)+2);
                    }
                    i+=2;
                }
                else{
                    v+=(codei-1);
                    i++;
                }
//                switch (get_code(i)) {
//                    case '\000':
//                        v--;
//                        i++;
//                        break;
//                    case '\001':
//                        i++;
//                        break;
//                    case '\002':
//                        v++;
//                        i++;
//                        break;
//                    case '\003':
//                        if (i == j - 1) {
//                            switch (get_code((i + 1))) {
//                                case '\000':
//                                    v += 2;
//                                    break;
//                                case '\001':
//                                    v += 3;
//                                    break;
//                                case '\002':
//                                    v += 4;
//                                    break;
//                                case '\003':
//                                    v += 5;
//                            }
//                        }
//                        i += 2;
//                }
            }
        }
        else{
            i=j;
            if (is_spike_code(i)) {//第一个delta可能是spikecode的后半部分。
                v+=(get_code(i)+2);
//                switch (get_code(i)) {
//                    case '\000':
//                        v += 2;
//                        break;
//                    case '\001':
//                        v += 3;
//                        break;
//                    case '\002':
//                        v += 4;
//                        break;
//                    case '\003':
//                        v += 5;
//                        break;
//                }
                i++;
            }
            while (i < N) {
                char codei= get_code(i);
                if(codei=='\003'){
                    if(i==N-1){
                        v-=(get_code(0)+2);
                    }
                    i+=2;
                }
                else{
                    v-=(codei-1);
                    i++;
                }
//                switch (get_code(i)) {
//                    case '\000':
//                        v++;
//                        i++;
//                        break;
//                    case '\001':
//                        i++;
//                        break;
//                    case '\002':
//                        v--;
//                        i++;
//                        break;
//                    case '\003':
//                        if (i == N - 1) {
//                            switch (get_code((0))) {
//                                case '\000':
//                                    v -= 2;
//                                    break;
//                                case '\001':
//                                    v -= 3;
//                                    break;
//                                case '\002':
//                                    v -= 4;
//                                    break;
//                                case '\003':
//                                    v -= 5;
//                            }
//                        }
//                        i += 2;
//                }
            }

        }
        return v;
    }
    bool get_extension(int j){
        return (((*(Se+j/8))>>(j%8))&1)=='\001';
//      return ((*(Se+j/8)>>(j%8))&1)=='\001';
//        int offset_byte=2*N+j;
//        return ((*(Sp+(offset_byte/8))>>(offset_byte%8))&1) == '\001';
    }
private:
    bool is_spike_code(int j){
        int count_11=0;
        for(int i=j-1;i>-1;i--){
            if(get_code(i)!='\003')
                return count_11 % 2 != 0;
            count_11++;
        }
        for(int i=N-1;i>j;i--){
            if(get_code(i)!='\003')
                return count_11 % 2 != 0;
            count_11++;
        }
        return true;
//        while(get_code((j-count_11-1<0)?(j-count_11-1+N):(j-count_11-1))=='\003'){
//            count_11++;
//        }
//        return count_11 % 2 != 0;
    }
//    int rho(unsigned __int128 x) const {
    static int rho(uint64_t x) {
        if(x<=3) return 31;
        uint32_t temp=(__builtin_clzll(x))/2;
        return  temp+ 1;
    }

    char get_code(int j){
        return  ((*(Sp+(j/4)))>>((2*j)%8))&3;
    }
    void set_code(int j,char code){
        if (code=='\000'){//00
            *(Sp+(j/4)) &= ~((unsigned char)(0b11<<(2*(j%4))));
        }
        else if(code=='\001'){//01
            *(Sp+(j/4)) |= (unsigned char)(0b1<<(2*(j%4)));
            *(Sp+(j/4)) &= ~((unsigned char)(0b1<<((2*(j%4))+1)));

        }
        else if(code=='\002'){//10
            *(Sp+(j/4)) &= ~((unsigned char)(0b1<<(2*(j%4))));
            *(Sp+(j/4)) |= (unsigned char)(0b1<<((2*(j%4))+1));
        }
        else{//11
            *(Sp+(j/4)) |= (unsigned char)(0b11<<(2*(j%4)));
        }
    }
    int get_delta(int j){
        int res;
        if(is_spike_code(j)) {
            res= -(get_code(j)+2);
//            switch (get_code(j)) {
//                case '\000':
//                    res = -2;
//                    break;
//                case '\001':
//                    res = -3;
//                    break;
//                case '\002':
//                    res = -4;
//                    break;
//                case '\003':
//                    res = -5;
//                    break;
//            }
        }
        else{
            char codej= get_code(j);
            if(codej=='\003'){
                res= get_code((j+1)%N)+2;
            }
            else{
                res=codej-1;
            }
//            switch (get_code(j)) {
//                case '\000':
//                    res=-1;
//                    break;
//                case '\001':
//                    res=0;
//                    break;
//                case '\002':
//                    res=1;
//                    break;
//                case '\003':
//                    switch (get_code((j+1)%N)){
//                        case '\000':
//                            res=2;
//                            break;
//                        case '\001':
//                            res=3;
//                            break;
//                        case '\002':
//                            res=4;
//                            break;
//                        case '\003':
//                            res=5;
//                    }
//            }
        }
        return res;
    }
    void set_extension(int j,int v){
        if(j>=N-1) return;
//        int offset_byte=2*N+j;
        if(v){
            *(Se+j/8) |= (unsigned char)(1<<(j%8));
        }
        else{
            *(Se+j/8) &= ~((unsigned char)(1<<(j%8)));
        }
    }
    bool tension(int j){
//        if(get_cell(j-1<0? j-1+N:j-1)<=-1 or get_cell(j)>=1){
//            return true;
//        }
//        return false;
        int jm1=j-1<0? j-1+N:j-1;
        if(is_spike_code(jm1) or  get_code(jm1)=='\000' or is_spike_code((j+1)%N) or ( !is_spike_code(j) and get_code(j)=='\002')){
            return true;
        }
        return false;
    }
    void up(int j, int up,int direction,int delta){
        //direction=0:left;1:right
        if(up<=0 && abs(delta)<2) return;

        int delta1;
        int k=0;
        if(direction==0){
            while(k<N-1){
                int jmk=(j-k<0)?N+j-k:j-k;
                int jmkm1=(j-k-1<0)?N+j-k-1:j-k-1;
                int  jmkm2=(j-k-2<0)?N+j-k-2:j-k-2;
                bool old_tension=this->tension(jmk);
                delta1= -delta;
                if(delta1>=2 && up-delta1<2){
                    delta=get_delta((j-k-3<0)?N+j-k-3:j-k-3);
                    if(j-k==0){
                        (*(Sbase))+=(up<<3);
                    }
                    if(up-delta1<=-2){
                        set_code(jmkm2,'\003');
                        switch (delta1-up) {
                            case 2:
                                set_code(jmkm1,'\000');
                                break;
                            case 3:
                                set_code(jmkm1,'\001');
                                break;
                            case 4:
                                set_code(jmkm1,'\002');
                                break;
                            case 5:
                                set_code(jmkm1,'\003');
                                break;
                        }
                    }
                    else{
                        set_code(jmkm2,'\002');
                        switch (delta1-up) {
                            case -1:
                                set_code(jmkm1,'\002');
                                up-=2;
                                break;
                            case 0:
                                set_code(jmkm1,'\001');
                                up-=1;
                                break;
                            case 1:
                                set_code(jmkm1,'\000');
                                break;
                        }
                    }
                    if (up<=0) return;
                    if(tension(jmk)){
                        set_extension(jmk,0);
                    }
                    else{
                        if(up==1 and (!old_tension or get_extension(jmk))){
                            set_extension(jmk,1);
                        }
                        else{
                            set_extension(jmk,0);
                        }
                    }
                    k+=2;

                }
                else{
                    if(j-k==0){
                        (*(Sbase))+=(up<<3);
                    }
                    delta=get_delta(jmkm2);
                    switch (up-delta1) {
                        case 0:
                            set_code(jmkm1,'\001');
                            break;
                        default:
                            set_code(jmkm1,'\002');
                    }
                    if (up!=0){
                        if(tension(jmk)){
                            set_extension(jmk,0);
                        }
                        else{
                            if(up==1 and (!old_tension or get_extension(jmk))){
                                set_extension(jmk,1);
                            }
                            else{
                                set_extension(jmk,0);
                            }
                        }
                    }
                    if (up-delta1-1>0){
                        up=up-delta1-1;
                        k+=1;
                    }
                    else{
                        break;
                    }
                }


            }
        }
        else{//向右传递。
            while(k<N-1){
                int jpk=(j+k)%N;
                int jpkp1=(j+k+1)%N;
                bool old_tension=this->tension(jpk);
                delta1= delta;//右侧第一个delta
                if(delta1>=2 && up-delta1<2){
                    /*下行为暂存右侧第三个delta。
                    在该分支中，我们可能会改变一个spikecode，将它的‘11’前缀替换成其他取值
                     我们在未改变其之前，将该delta解码出来保存，并在下次循环时将它赋值给delta1*/
                    delta=get_delta((j+k+2)%N);
                    if(j+k==N){
                        (*(Sbase))+=(up<<3);//如果此时要处理的cell为第一个cell，我们不仅要更改delta0(绝对index)和delta40，还要同步更新base
                    }

                    if(up-delta1<=-2){//该分支表示经过up后，delta1和delta2仍是一个spikecode的情况。
                        set_code(jpk,'\003');//设置前缀
//                        set_code(jpkp1,(char)(delta1-up-2));
                        switch (up-delta1) {//设置spikecode的后半部分
                            case -2:
                                set_code(jpkp1,'\000');
                                break;
                            case -3:
                                set_code(jpkp1,'\001');
                                break;
                            case -4:
                                set_code(jpkp1,'\002');
                                break;
                            case -5:
                                set_code(jpkp1,'\003');
                                break;
                        }
                    }
                    else{//经过up后，delta1和delta2不再是是一个spikecode。
                        set_code(jpkp1,'\000');//delta1一定是-1
//                        set_code(jpk,(char)(delta1-up+1));
//                        up=delta1-1;
                        switch (up-delta1) {
                            case -1:
                                set_code(jpk,'\002');
                                break;
                            case 0:
                                set_code(jpk,'\001');
                                up-=1;
                                break;
                            case 1:
                                set_code(jpk,'\000');
                                up-=2;
                                break;
                        }
                    }
                    if (up<=0) return;
                    if(tension(jpk)){
                        set_extension(jpk,0);
                    }
                    else{
                        if(up==1 and (!old_tension or get_extension(jpk))){
                            set_extension(jpk,1);
                        }
                        else{
                            set_extension(jpk,0);
                        }
                    }

                    k+=2;

                }
                else{
                    delta=get_delta(jpkp1);
                    if(j+k==N ){
                        (*(Sbase))+=(up<<3);//如果此时要处理的cell为第一个cell，我们不仅要更改delta0(绝对index)和delta40，还要同步更新base
                    }
//                    set_code(jpk,(up-delta1)?'\000':'\001');
                    switch (up-delta1) {
                        case 0:
                            set_code(jpk,'\001');
                            break;
                        default:
                            set_code(jpk,'\000');
                    }
                    if(up!=0){
                        if(tension(jpk)){
                            set_extension(jpk,0);
                        }
                        else{
                            if(up==1 and (!old_tension or get_extension(jpk))){
                                set_extension(jpk,1);
                            }
                            else{
                                set_extension(jpk,0);
                            }
                        }
                    }
                    if (up-delta1-1>0){
                        up=up-delta1-1;
                        k+=1;
                    }
                    else{
                        break;
                    }
                }
            }
        }
    }
    void adjust(int j,int up){
        int new_delta1,new_delta2,min_delta,jm1;
        jm1=(j-1<0)?N+j-1:j-1;
        new_delta1=up+get_delta(jm1);
        new_delta2=up-get_delta(j);
        min_delta=min(new_delta1,new_delta2);

        bool old_tension=this->tension(j);
        bool old_extension=this->get_extension(j);

        if(min_delta>VLIMIT){
            int delta0,delta3,jp1,actual_delta;
            jp1=(j+1)%N;
            actual_delta=min_delta<SLIMIT?min_delta:SLIMIT;
            delta0=get_delta((j-2<0)?N+j-2:j-2);
            delta3=get_delta(jp1);
            set_code(jm1,'\003');
            set_code(j,(char)(min_delta<5?min_delta-2:3));
//            switch (min_delta) {
//                case 2:
//                    set_code(j,'\000');
//                    break;
//                case 3:
//                    set_code(j,'\001');
//                    break;
//                case 4:
//                    set_code(j,'\002');
//                    break;
//                case 5:
//                    set_code(j,'\003');
//                    break;
//                default:
//                    set_code(j,'\003');
//            }
            if(j==0){
                (*(Sbase))+=(up<<3);
            }
            if(j==N-1){
                (*(Sbase))+=((new_delta2-actual_delta)<<3);
            }
            this->up(jm1,new_delta1-actual_delta,0,delta0);
            this->up(jp1,new_delta2-actual_delta,1,delta3);
        }
        else{
            this->up(j,up,0,new_delta1-up);
            this->up(j,up,1,-new_delta2+up);
        }

        if(tension(j)){
            set_extension(j,1);
        }
        else{
            if(up==1 and (!old_tension or old_extension)){
                set_extension(j,1);
            }
            else{
                set_extension(j,0);
            }
        }
    }
};
