/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

#ifndef CPC_SKETCH_IMPL_HPP_
#define CPC_SKETCH_IMPL_HPP_

#include <stdexcept>
#include <cmath>
#include <cstring>
#include <sstream>

#include "cpc_confidence.hpp"
#include "kxp_byte_lookup.hpp"
#include "inv_pow2_table.hpp"
#include "cpc_util.hpp"
#include "icon_estimator.hpp"
#include "serde.hpp"
#include "count_zeros.hpp"

namespace datasketches {

template<typename A>
void cpc_init() {
  get_compressor<A>(); // this initializes a global static instance of the compressor on the first use
}

template<typename A>
cpc_sketch_alloc<A>::cpc_sketch_alloc(uint32_t k, uint64_t seed, const A& allocator):
k(k),
seed(seed),
was_merged(true),
num_coupons(0),
surprising_value_table(2, k<2? 6+1:6 + ceil(std::log2((double)k)), allocator),
sliding_window(allocator),
window_offset(0),
first_interesting_column(0),
kxp(k),
hip_est_accum(0)
{
  check_k(k);
}

template<typename A>
A cpc_sketch_alloc<A>::get_allocator() const {
  return sliding_window.get_allocator();
}

template<typename A>
uint8_t cpc_sketch_alloc<A>::get_k() const {
  return k;
}

template<typename A>
bool cpc_sketch_alloc<A>::is_empty() const {
  return num_coupons == 0;
}

template<typename A>
double cpc_sketch_alloc<A>::get_estimate() const {
//  if (!was_merged) return get_hip_estimate();
  return get_icon_estimate();
}

template<typename A>
double cpc_sketch_alloc<A>::get_hip_estimate() const {
  return hip_est_accum;
}

template<typename A>
double cpc_sketch_alloc<A>::get_icon_estimate() const {
  return compute_icon_estimate(k, num_coupons);
}

//template<typename A>
//double cpc_sketch_alloc<A>::get_lower_bound(unsigned kappa) const {
//  if (kappa < 1 || kappa > 3) {
//    throw std::invalid_argument("kappa must be 1, 2 or 3");
//  }
//  if (!was_merged) return get_hip_confidence_lb<A>(*this, kappa);
//  return get_icon_confidence_lb<A>(*this, kappa);
//}
//
//template<typename A>
//double cpc_sketch_alloc<A>::get_upper_bound(unsigned kappa) const {
//  if (kappa < 1 || kappa > 3) {
//    throw std::invalid_argument("kappa must be 1, 2 or 3");
//  }
//  if (!was_merged) return get_hip_confidence_ub<A>(*this, kappa);
//  return get_icon_confidence_ub<A>(*this, kappa);
//}

//template<typename A>
//void cpc_sketch_alloc<A>::update(const std::string& value) {
//  if (value.empty()) return;
//  update(value.c_str(), value.length());
//}

template<typename A>
uint8_t cpc_sketch_alloc<A>::update(uint64_t value) {
  update(&value, sizeof(value));

  return surprising_value_table.get_lg_size();
}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(int64_t value) {
//  update(&value, sizeof(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(uint32_t value) {
//  update(static_cast<int32_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(int32_t value) {
//  update(static_cast<int64_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(uint16_t value) {
//  update(static_cast<int16_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(int16_t value) {
//  update(static_cast<int64_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(uint8_t value) {
//  update(static_cast<int8_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(int8_t value) {
//  update(static_cast<int64_t>(value));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(double value) {
//  union {
//    int64_t long_value;
//    double double_value;
//  } ldu;
//  if (value == 0.0) {
//    ldu.double_value = 0.0; // canonicalize -0.0 to 0.0
//  } else if (std::isnan(value)) {
//    ldu.long_value = 0x7ff8000000000000L; // canonicalize NaN using value from Java's Double.doubleToLongBits()
//  } else {
//    ldu.double_value = value;
//  }
//  update(&ldu, sizeof(ldu));
//}
//
//template<typename A>
//void cpc_sketch_alloc<A>::update(float value) {
//  update(static_cast<double>(value));
//}

static inline uint32_t row_col_from_two_hashes(uint64_t hash0, uint64_t hash1, uint32_t k) {
  if (k > (1<<26)) throw std::logic_error("lg_k > 26");
//  const uint32_t k = 1 << lg_k;
  uint8_t col = count_leading_zeros_in_u64(hash1); // 0 <= col <= 64
  if (col > 63) col = 63; // clip so that 0 <= col <= 63
  const uint32_t row = hash0 % k;
  uint32_t row_col = (row << 6) | col;
  // To avoid the hash table's "empty" value, we change the row of the following pair.
  // This case is extremely unlikely, but we might as well handle it.
  if (row_col == UINT32_MAX) row_col ^= 1 << 6;
  return row_col;
}

template<typename A>
void cpc_sketch_alloc<A>::update(const void* value, size_t size) {
  HashState hashes;
  MurmurHash3_x64_128(value, size, seed, hashes);
  row_col_update(row_col_from_two_hashes(hashes.h1, hashes.h2, k));
}

template<typename A>
void cpc_sketch_alloc<A>::row_col_update(uint32_t row_col) {
  const uint8_t col = row_col & 63;
  if (col < first_interesting_column) return; // important speed optimization
  // window size is 0 until sketch is promoted from sparse to windowed
  if (sliding_window.size() == 0) {
    update_sparse(row_col);
  } else {
    update_windowed(row_col);
  }
}

template<typename A>
void cpc_sketch_alloc<A>::update_sparse(uint32_t row_col) {
//  const uint32_t k = 1 << lg_k;
  const uint64_t c32pre = static_cast<uint64_t>(num_coupons) << 5;
  if (c32pre >= 3 * k) throw std::logic_error("c32pre >= 3 * k"); // C < 3K/32, in other words flavor == SPARSE
  bool is_novel = surprising_value_table.maybe_insert(row_col);
  if (is_novel) {
    num_coupons++;
    update_hip(row_col);
    const uint64_t c32post = static_cast<uint64_t>(num_coupons) << 5;
    if (c32post >= 3 * k) promote_sparse_to_windowed(); // C >= 3K/32
  }
}

// the flavor is HYBRID, PINNED, or SLIDING
template<typename A>
void cpc_sketch_alloc<A>::update_windowed(uint32_t row_col) {
  if (window_offset > 56) throw std::logic_error("wrong window offset");
//  const uint32_t k = 1 << lg_k;
  const uint64_t c32pre = static_cast<uint64_t>(num_coupons) << 5;
  if (c32pre < 3 * k) throw std::logic_error("c32pre < 3 * k"); // C < 3K/32, in other words flavor >= HYBRID
  const uint64_t c8pre = static_cast<uint64_t>(num_coupons) << 3;
  const uint64_t w8pre = static_cast<uint64_t>(window_offset) << 3;
  if (c8pre >= (27 + w8pre) * k) throw std::logic_error("c8pre is wrong"); // C < (K * 27/8) + (K * window_offset)

  bool is_novel = false;
  const uint8_t col = row_col & 63;

  if (col < window_offset) { // track the surprising 0's "before" the window
    is_novel = surprising_value_table.maybe_delete(row_col); // inverted logic
  } else if (col < window_offset + 8) { // track the 8 bits inside the window
    if (col < window_offset) throw std::logic_error("col < window_offset");
    const uint32_t row = row_col >> 6;
    const uint8_t old_bits = sliding_window[row];
    const uint8_t new_bits = old_bits | (1 << (col - window_offset));
    if (new_bits != old_bits) {
      sliding_window[row] = new_bits;
      is_novel = true;
    }
  } else { // track the surprising 1's "after" the window
    if (col < window_offset + 8) throw std::logic_error("col < window_offset + 8");
    is_novel = surprising_value_table.maybe_insert(row_col); // normal logic
  }

  if (is_novel) {
    num_coupons++;
    update_hip(row_col);
    const uint64_t c8post = static_cast<uint64_t>(num_coupons) << 3;
    if (c8post >= (27 + w8pre) * k) {
      move_window();
      if (window_offset < 1 || window_offset > 56) throw std::logic_error("wrong window offset");
      const uint64_t w8post = static_cast<uint64_t>(window_offset) << 3;
      if (c8post >= (27 + w8post) * k) throw std::logic_error("c8pre is wrong"); // C < (K * 27/8) + (K * window_offset)
    }
  }
}

// Call this whenever a new coupon has been collected.
template<typename A>
void cpc_sketch_alloc<A>::update_hip(uint32_t row_col) {
//  const uint32_t k = 1 << lg_k;
  const uint8_t col = row_col & 63;
  const double one_over_p = static_cast<double>(k) / kxp;
  hip_est_accum += one_over_p;
  kxp -= INVERSE_POWERS_OF_2[col + 1]; // notice the "+1"
}

// In terms of flavor, this promotes SPARSE to HYBRID
template<typename A>
void cpc_sketch_alloc<A>::promote_sparse_to_windowed() {
//  const uint32_t k = 1 << lg_k;
  const uint64_t c32 = static_cast<uint64_t>(num_coupons) << 5;
//  if (!(c32 == 3 * k || (lg_k == 4 && c32 > 3 * k))) throw std::logic_error("wrong c32");

  sliding_window.resize(k, 0); // zero the memory (because we will be OR'ing into it)

  u32_table<A> new_table(2, k<2? 6+1:6 + ceil(std::log2((double)k)), sliding_window.get_allocator());

  const uint32_t* old_slots = surprising_value_table.get_slots();
  const uint32_t old_num_slots = 1 << surprising_value_table.get_lg_size();

  if (window_offset != 0) throw std::logic_error("window_offset != 0");

  for (uint32_t i = 0; i < old_num_slots; i++) {
    const uint32_t row_col = old_slots[i];
    if (row_col != UINT32_MAX) {
      const uint8_t col = row_col & 63;
      if (col < 8) {
        const uint32_t row = row_col >> 6;
        sliding_window[row] |= 1 << col;
      } else {
        // cannot use u32_table::must_insert(), because it doesn't provide for growth
        const bool is_novel = new_table.maybe_insert(row_col);
        if (!is_novel) throw std::logic_error("is_novel != true");
      }
    }
  }

  surprising_value_table = std::move(new_table);
}

template<typename A>
void cpc_sketch_alloc<A>::move_window() {
  const uint8_t new_offset = window_offset + 1;
  if (new_offset > 56) throw std::logic_error("new_offset > 56");
//  if (new_offset != determine_correct_offset(lg_k, num_coupons)) throw std::logic_error("new_offset is wrong");

  if (sliding_window.size() == 0) throw std::logic_error("no sliding window");
//  const uint32_t k = 1 << lg_k;

  // Construct the full-sized bit matrix that corresponds to the sketch
  vector_u64<A> bit_matrix = build_bit_matrix();

  // refresh the KXP register on every 8th window shift.
  if ((new_offset & 0x7) == 0) refresh_kxp(bit_matrix.data());

  surprising_value_table.clear(); // the new number of surprises will be about the same

  const uint64_t mask_for_clearing_window = (static_cast<uint64_t>(0xff) << new_offset) ^ UINT64_MAX;
  const uint64_t mask_for_flipping_early_zone = (static_cast<uint64_t>(1) << new_offset) - 1;
  uint64_t all_surprises_ored = 0;

  for (uint32_t i = 0; i < k; i++) {
    uint64_t pattern = bit_matrix[i];
    sliding_window[i] = (pattern >> new_offset) & 0xff;
    pattern &= mask_for_clearing_window;
    // The following line converts surprising 0's to 1's in the "early zone",
    // (and vice versa, which is essential for this procedure's O(k) time cost).
    pattern ^= mask_for_flipping_early_zone;
    all_surprises_ored |= pattern; // a cheap way to recalculate first_interesting_column
    while (pattern != 0) {
      const uint8_t col = count_trailing_zeros_in_u64(pattern);
      pattern = pattern ^ (static_cast<uint64_t>(1) << col); // erase the 1
      const uint32_t row_col = (i << 6) | col;
      const bool is_novel = surprising_value_table.maybe_insert(row_col);
      if (!is_novel) throw std::logic_error("is_novel != true");
    }
  }

  window_offset = new_offset;

  first_interesting_column = count_trailing_zeros_in_u64(all_surprises_ored);
  if (first_interesting_column > new_offset) first_interesting_column = new_offset; // corner case
}

// The KXP register is a double with roughly 50 bits of precision, but
// it might need roughly 90 bits to track the value with perfect accuracy.
// Therefore we recalculate KXP occasionally from the sketch's full bitmatrix
// so that it will reflect changes that were previously outside the mantissa.
template<typename A>
void cpc_sketch_alloc<A>::refresh_kxp(const uint64_t* bit_matrix) {
//  const uint32_t k = 1 << lg_k;

  // for improved numerical accuracy, we separately sum the bytes of the U64's
  double byte_sums[8]; // allocating on the stack
  std::fill(byte_sums, &byte_sums[8], 0);

  for (size_t i = 0; i < k; i++) {
    uint64_t word = bit_matrix[i];
    for (unsigned j = 0; j < 8; j++) {
      const uint8_t byte = word & 0xff;
      byte_sums[j] += KXP_BYTE_TABLE[byte];
      word >>= 8;
    }
  }

  double total = 0.0;
  for (int j = 7; j >= 0; j--) { // the reverse order is important
    const double factor = INVERSE_POWERS_OF_2[8 * j]; // pow (256.0, (-1.0 * ((double) j)));
    total += factor * byte_sums[j];
  }

  kxp = total;
}

//template<typename A>
//string<A> cpc_sketch_alloc<A>::to_string() const {
//  // Using a temporary stream for implementation here does not comply with AllocatorAwareContainer requirements.
//  // The stream does not support passing an allocator instance, and alternatives are complicated.
//  std::ostringstream os;
//  os << "### CPC sketch summary:" << std::endl;
//  os << "   k           : " << std::to_string(k) << std::endl;
//  os << "   seed hash      : " << std::hex << compute_seed_hash(seed) << std::dec << std::endl;
//  os << "   C              : " << num_coupons << std::endl;
//  os << "   flavor         : " << determine_flavor() << std::endl;
//  os << "   merged         : " << (was_merged ? "true" : "false") << std::endl;
//  if (!was_merged) {
//    os << "   HIP estimate   : " << hip_est_accum << std::endl;
//    os << "   kxp            : " << kxp << std::endl;
//  }
//  os << "   interesting col: " << std::to_string(first_interesting_column) << std::endl;
//  os << "   table entries  : " << surprising_value_table.get_num_items() << std::endl;
//  os << "   window         : " << (sliding_window.size() == 0 ? "not " : "") <<  "allocated" << std::endl;
//  if (sliding_window.size() > 0) {
//    os << "   window offset  : " << std::to_string(window_offset) << std::endl;
//  }
//  os << "### End sketch summary" << std::endl;
//  return string<A>(os.str().c_str(), sliding_window.get_allocator());
//}


///*
// * These empirical values for the 99.9th percentile of size in bytes were measured using 100,000
// * trials. The value for each trial is the maximum of 5*16=80 measurements that were equally
// * spaced over values of the quantity C/K between 3.0 and 8.0. This table does not include the
// * worst-case space for the preamble, which is added by the function.
// */
//static const uint8_t CPC_EMPIRICAL_SIZE_MAX_LGK = 19;
//static const size_t CPC_EMPIRICAL_MAX_SIZE_BYTES[]  = {
//    24,     // lg_k = 4
//    36,     // lg_k = 5
//    56,     // lg_k = 6
//    100,    // lg_k = 7
//    180,    // lg_k = 8
//    344,    // lg_k = 9
//    660,    // lg_k = 10
//    1292,   // lg_k = 11
//    2540,   // lg_k = 12
//    5020,   // lg_k = 13
//    9968,   // lg_k = 14
//    19836,  // lg_k = 15
//    39532,  // lg_k = 16
//    78880,  // lg_k = 17
//    157516, // lg_k = 18
//    314656  // lg_k = 19
//};
//static const double CPC_EMPIRICAL_MAX_SIZE_FACTOR = 0.6; // 0.6 = 4.8 / 8.0
//static const size_t CPC_MAX_PREAMBLE_SIZE_BYTES = 40;

//template<typename A>
//size_t cpc_sketch_alloc<A>::get_max_serialized_size_bytes(uint8_t lg_k) {
//  check_lg_k(lg_k);
//  if (lg_k <= CPC_EMPIRICAL_SIZE_MAX_LGK) return CPC_EMPIRICAL_MAX_SIZE_BYTES[lg_k - CPC_MIN_LG_K] + CPC_MAX_PREAMBLE_SIZE_BYTES;
//  const uint32_t k = 1 << lg_k;
//  return (int) (CPC_EMPIRICAL_MAX_SIZE_FACTOR * k) + CPC_MAX_PREAMBLE_SIZE_BYTES;
//}

template<typename A>
void cpc_sketch_alloc<A>::check_k(uint32_t k) {
//  if (k < (1<<CPC_MIN_LG_K) || k > (1<<CPC_MAX_LG_K)) {
//    throw std::invalid_argument("lg_k must be >= " + std::to_string(CPC_MIN_LG_K) + " and <= " + std::to_string(CPC_MAX_LG_K) + ": " + std::to_string(k));
//  }
}

template<typename A>
uint32_t cpc_sketch_alloc<A>::get_num_coupons() const {
  return num_coupons;
}

//template<typename A>
//bool cpc_sketch_alloc<A>::validate() const {
//  vector_u64<A> bit_matrix = build_bit_matrix();
//  const uint64_t num_bits_set = count_bits_set_in_matrix(bit_matrix.data(), 1ULL << lg_k);
//  return num_bits_set == num_coupons;
//}

//template<typename A>
//cpc_sketch_alloc<A>::cpc_sketch_alloc(uint8_t lg_k, uint32_t num_coupons, uint8_t first_interesting_column,
//    u32_table<A>&& table, vector_u8<A>&& window, bool has_hip, double kxp, double hip_est_accum, uint64_t seed):
//lg_k(lg_k),
//seed(seed),
//was_merged(!has_hip),
//num_coupons(num_coupons),
//surprising_value_table(std::move(table)),
//sliding_window(std::move(window)),
//window_offset(determine_correct_offset(lg_k, num_coupons)),
//first_interesting_column(first_interesting_column),
//kxp(kxp),
//hip_est_accum(hip_est_accum)
//{}
//
//template<typename A>
//uint8_t cpc_sketch_alloc<A>::get_preamble_ints(uint32_t num_coupons, bool has_hip, bool has_table, bool has_window) {
//  uint8_t preamble_ints = 2;
//  if (num_coupons > 0) {
//    preamble_ints += 1; // number of coupons
//    if (has_hip) {
//      preamble_ints += 4; // HIP
//    }
//    if (has_table) {
//      preamble_ints += 1; // table data length
//      // number of values (if there is no window it is the same as number of coupons)
//      if (has_window) {
//        preamble_ints += 1;
//      }
//    }
//    if (has_window) {
//      preamble_ints += 1; // window length
//    }
//  }
//  return preamble_ints;
//}

template<typename A>
typename cpc_sketch_alloc<A>::flavor cpc_sketch_alloc<A>::determine_flavor() const {
  return determine_flavor(k, num_coupons);
}

template<typename A>
typename cpc_sketch_alloc<A>::flavor cpc_sketch_alloc<A>::determine_flavor(uint32_t k, uint64_t c) {
  const uint64_t c2 = c << 1;
  const uint64_t c8 = c << 3;
  const uint64_t c32 = c << 5;
  if (c == 0)      return EMPTY;    //    0  == C <    1
  if (c32 < 3 * k) return SPARSE;   //    1  <= C <   3K/32
  if (c2 < k)      return HYBRID;   // 3K/32 <= C <   K/2
  if (c8 < 27 * k) return PINNED;   //   K/2 <= C < 27K/8
  else             return SLIDING;  // 27K/8 <= C
}

//template<typename A>
//uint8_t cpc_sketch_alloc<A>::determine_correct_offset(uint8_t lg_k, uint64_t c) {
//  const uint32_t k = 1 << lg_k;
//  const int64_t tmp = static_cast<int64_t>(c << 3) - static_cast<int64_t>(19 * k); // 8C - 19K
//  if (tmp < 0) return 0;
//  return static_cast<uint8_t>(tmp >> (lg_k + 3)); // tmp / 8K
//}

template<typename A>
vector_u64<A> cpc_sketch_alloc<A>::build_bit_matrix() const {
//  const uint32_t k = this->k;
  if (window_offset > 56) throw std::logic_error("offset > 56");

  // Fill the matrix with default rows in which the "early zone" is filled with ones.
  // This is essential for the routine's O(k) time cost (as opposed to O(C)).
  const uint64_t default_row = (static_cast<uint64_t>(1) << window_offset) - 1;
  vector_u64<A> matrix(k, default_row, sliding_window.get_allocator());

  if (num_coupons == 0) return matrix;

  if (sliding_window.size() > 0) { // In other words, we are in window mode, not sparse mode
    for (size_t i = 0; i < k; i++) { // set the window bits, trusting the sketch's current offset
      matrix[i] |= static_cast<uint64_t>(sliding_window[i]) << window_offset;
    }
  }

  const uint32_t* slots = surprising_value_table.get_slots();
  const uint32_t num_slots = 1 << surprising_value_table.get_lg_size();
  for (size_t i = 0; i < num_slots; i++) {
    const uint32_t row_col = slots[i];
    if (row_col != UINT32_MAX) {
      const uint8_t col = row_col & 63;
      const uint32_t row = row_col >> 6;
      // Flip the specified matrix bit from its default value.
      // In the "early" zone the bit changes from 1 to 0.
      // In the "late" zone the bit changes from 0 to 1.
      matrix[row] ^= static_cast<uint64_t>(1) << col;
    }
  }
  return matrix;
}

} /* namespace datasketches */

#endif
