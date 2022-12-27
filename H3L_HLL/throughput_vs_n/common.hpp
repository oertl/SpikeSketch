#ifndef HYPERLOGLOGLOG_COMMON
#define HYPERLOGLOGLOG_COMMON

//#include <arpa/inet.h>
#include <cstdint>

namespace hyperlogloglog {

    template<typename T>
    inline int clz(T x);

    template<>
    inline int clz(unsigned int x) {
        return __builtin_clz(x);
    }

    template<>
    inline int clz(unsigned long x) {
        return __builtin_clzl(x);
    }

    template<>
    inline int clz(unsigned long long x) {
        return __builtin_clzll(x);
    }

    template<typename T>
    int rho(T x) {
        if(x<=1) return 63;
        return clz(x) + 1;
    }

//
//    int rho(uint64_t x) {
//        for (int i=0; i < 63; i++){
//            if(x&1==1){
//                return i+1;
//            }
//            x=x>>1;
//        }
//        return 63;
//    }

    template<typename T>
    constexpr T log2i(T x) {
        return x<2 ? 1 : (int)ceil(log2((double)x));
    }
}

#endif // HYPERLOGLOGLOG_COMMON
