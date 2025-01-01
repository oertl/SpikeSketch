#include <iostream>
#include <random>
#include <vector>

#include "SpikeSketch/spike_sketch_extend.h"

void merge(vector<vector<spike_sketch>> spikeSketchArray,
           vector<spike_sketch> &merged_ss);

double mutiBktQuery(vector<spike_sketch> &spikeSketchArray, double alpha0,
                    double alpha1, double beta0, double beta1, double coe);

const int n = 20;    // Number of cells in a single bucket
const int ncode = 4; // Number of bits in a single cell
const constexpr int p = 12;
const uint32_t numOfBuckets = 128;
const uint32_t seed = 0x529b9601;

vector<spike_sketch> create() {
  vector<spike_sketch> spikeSketchArray;
  spikeSketchArray.reserve(numOfBuckets);
  for (uint32_t bktIdx = 0; bktIdx < numOfBuckets; bktIdx++) {
    spikeSketchArray.emplace_back(n, p, ncode, seed);
  }
  return spikeSketchArray;
}

void add(vector<spike_sketch> &sketch, uint64_t hash) {
  uint32_t tempInt32 = 0;
  MurmurHash3_x86_32(&hash, 8, seed + 231321, &tempInt32);
  spike_sketch ss = sketch[tempInt32 % numOfBuckets];
  ss.update(hash);
}

double estimate(const vector<spike_sketch> &sketch) {

  double alpha0 = 0.1;
  double alpha1 = 0.88;
  double beta0 = 1.12;
  double beta1 = 1.46;
  double myCoe = 0.573; // Coefficient of correction

  double estimate = mutiBktQuery(const_cast<vector<spike_sketch> &>(sketch),
                                 alpha0, alpha1, beta0, beta1, myCoe);
  return estimate;
}

int main() {

  mt19937_64 rng(4);

  uint64_t numElements = 1000;

  vector<spike_sketch> sketch1 = create();
  for (uint64_t c = 0; c < numElements; ++c) {
    add(sketch1, rng());
  }
  vector<spike_sketch> sketch2 = create();
  for (uint64_t c = 0; c < numElements; ++c) {
    add(sketch2, rng());
  }
  vector<spike_sketch> merged_sketch = create();
  vector<vector<spike_sketch>> sketches_to_merge = {sketch1, sketch2};
  merge(sketches_to_merge, merged_sketch);

  cout << "distinct count 1 = " << estimate(sketch1) << endl;
  cout << "distinct count 2 = " << estimate(sketch2) << endl;
  cout << "distinct count merged = " << estimate(merged_sketch) << endl;
}