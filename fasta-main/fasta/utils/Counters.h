/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <chrono>
#include <cstdint>
#include <sstream>
#include <fmt/format.h>

namespace counter {

  /**
   * Counter class
   */
  typedef struct Counter {
    uint_fast64_t count = 0;
    bool enabled = true;

    inline void reset(){
      count = 0;
    }
    Counter& operator+= (const Counter& rhs){
      count += rhs.count;
      return *this;
    }
    Counter& operator+= (const uint_fast64_t rhs){
      count += rhs;
      return *this;
    }
    Counter& operator++(){
      count++;
      return *this;
    }
    Counter& operator++(int){
      auto temp = *this;
      ++*this;
      return temp;
    }
    Counter& operator=(const uint64_t rhs){
      count = rhs;
      return *this;
    }
    friend Counter operator+ (Counter lhs, const Counter& rhs){
      lhs.count += rhs.count;
      return lhs;
    }
  } Counter;


  std::ostream& operator<< (std::ostream &os, const Counter &c) {
    return (os << c.count);
  }

  /** Counter types (index in Counters array) */
  enum CounterTypes {
    INDEX_SIZE,     // number of non-zeros in index
    NUM_CANDIDATES,   // number of candidates
    PRUNE_CG_L2,    // number of candidates pruned during candidate generation
    PRUNE_CV_L2,    // number of candidates pruned during candidate verification
    DOT_PRODUCTS,   // number of full dot products executed
    RESULTS,      // number of results
    NCOUNTERS
  };

  /**
   * Print out array of counters
   * @param  counters Array of counters to print
   * @param  prefix   Prefix for each timer line.
   * @return      Human readable counters.
   */
  std::string get_counts_string(const Counter* counters, std::string prefix=""){
    std::ostringstream s;

    s << prefix << "Index size: " << counters[INDEX_SIZE] << std::endl;
    s << prefix << "Number of candidates: " << counters[NUM_CANDIDATES] << std::endl;
    #ifdef EXTRA_COUNTERS
    if(counters[PRUNE_CG_L2].count){
      s << prefix << "Candidate Generation L2 prunning: " << counters[PRUNE_CG_L2] << std::endl;
    }
    if(counters[PRUNE_CV_L2].count){
      s << prefix << "Candidate Verification L2 prunning: " << counters[PRUNE_CV_L2] << std::endl;
    }
    #endif
    s << prefix << "Number of full dot-products: " << counters[DOT_PRODUCTS] << std::endl;
    s << prefix << "Result size: " << counters[RESULTS] << std::endl;
    return s.str();
  }

  /**
   * Add up the counters from source to target.
   * @param source Source timers array.
   * @param target Target timers array.
   */
  void add_counters(
    const Counter* const source,
    Counter* const target
  ){
    for(uint8_t i=0; i < NCOUNTERS; ++i){
      target[i].count += source[i].count;
    }
  }

}  // namespace counter
