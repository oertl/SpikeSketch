# SpikeSketch
SpikeSketch is a novel cardinality estimator.It will be presented at the IEEE International Conference on Computer Communications ([INFOCOM 2023](https://infocom2023.ieee-infocom.org/)). This is the repository for SpikeSketch's source codes.
## Introduction
Cardinality estimation is a fundamental problem with diverse practical applications. HyperLogLog (HLL) has become a standard in practice because it offers good memory efficiency, constant update time, and mergeability. Some recent work achieved better memory efficiency, but typically at the cost of impractical update time or losing mergeability, making them incompatible with applications like network-wide traffic measurement. This work presents SpikeSketch, a better cardinality estimator that reduces memory usage of HLL by 37% without sacrificing other crucial metrics. We adopt a bucket-based data structure to promise constant update time, design a smoothed log$_4$ ranking and a spike coding scheme to compress cardinality observables into buckets, and propose a lightweight mergeable lossy compression to balance memory usage, information loss, and mergeability. Then we derive an unbiased estimator for recovering cardinality from the lossy-compressed sketch. Theoretical and empirical results show that SpikeSketch can work as a drop-in replacement for HLL because it achieves a near-optimal MVP (memory-variance-product) of 4.08 (37% smaller than HLL) with constant update time and mergeability. Its memory efficiency even defeats ACPC and HLLL, the state-of-the-art lossless-compressed sketches using linear-time compression to reduce memory usage. 
## About this repo
In our paper, we evaluated the performance of our proposed SpikeSketch and compared it with state-of-the-art algorithms, including HyperLogLog (HLL) , ACPC, and HyperLogLogLog (HLLL). The implementation code of these four algorithms are given in this repository.
## How to run these codes?
We strongly recommend that you compile these code using CMake. We have given the corresponding CMakeLists.txt in each directory. The compiling command is as follows.
```shell
cmake .
cmake --build . --config Release
```
