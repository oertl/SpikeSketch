g++ -Wl,--allow-multiple-definition test_merge.cpp SpikeSketch/main_extend.cpp SpikeSketch/utils/MurmurHash3.cpp SpikeSketch/impl/ss_query1.cpp SpikeSketch/merge.cpp -o testMerge
./testMerge