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
#include <iostream>
#include <sstream>
#include <fmt/format.h>

namespace timer {

  /**
   * Timer class
   */
  typedef struct Timer {
    uint64_t time = 0;
    bool started = false;
    std::chrono::time_point<std::chrono::high_resolution_clock> timer;
    inline bool active() const{
      return time || started;
    }
    inline void start(){
      timer = std::chrono::high_resolution_clock::now();
      started = true;
    }
    inline void stop(){
      started = false;
      auto end = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - timer);
      time += elapsed.count();
    }
    inline uint64_t get_time(){
      if(started){
        stop();
      }
      return time;
    }
    inline double get_time_seconds(){
      if(started){
        stop();
      }
      return (double) ((long double) time * 1e-9);
    }
    inline void reset(){
      time = 0;
    }
    inline std::string get_time_string() const{
      double seconds;
      uint32_t minutes, hours, days, years;
      seconds = (double) ((long double) time * 1e-9);
      if(seconds < 60.0){
        return fmt::format("{:0.4f}", seconds);
      }
      minutes = (uint32_t) (seconds / 60.0);
      seconds -= minutes * 60;
      if(minutes < 60){
        return fmt::format("{:02d}:{:0.4f}", minutes, seconds);
      }
      hours = minutes / 60;
      minutes -= hours * 60;
      if(hours < 24){
        return fmt::format("{:02d}:{:02d}:{:0.4f}", hours, minutes, seconds);
      }
      days = hours / 24;
      hours -= days * 24;
      if(days < 365){
        return fmt::format("{:03d}:{:02d}:{:02d}:{:0.4f}", days, hours, minutes, seconds);
      }
      years = days / 365;
      days -= years * 365;
      return fmt::format("{}:{:03d}:{:02d}:{:02d}:{:0.4f}", years, days, hours, minutes, seconds);
    }
    Timer& operator+= (const Timer& rhs){
      if(started){
        stop();
      }
      time += rhs.time;
      return *this;
    }
    Timer& operator/= (const uint64_t rhs){
      if(started){
        stop();
      }
      time /= rhs;
      return *this;
    }
    friend Timer operator+ (Timer lhs, const Timer& rhs){
      if(lhs.started){
        lhs.stop();
      }
      lhs.time += rhs.time;
      return lhs;
    }
  } Timer;

  std::ostream& operator<< (std::ostream &os, const Timer &t) {
    return (os << t.get_time_string());
  }

  /** Timer types (index in Timers array) */
  enum TimerTypes {
    EXEC,
    READ,
    WRITE,
    PREPROCESS,
    QUERY_PREPROCESS,
    SEARCH,
    CANDIDATE_GENERATION,
    CANDIDATE_VERIFICATION,
    NTIMERS
  };

  /**
   * Create a print-out string for the set of timers.
   * @param  timers Array of timers
   * @param  prefix Prefix for each timer line.
   * @return    Human readable timers.
   */
  std::string get_times_string(
    const Timer* timers,
    std::string prefix=""
  ){
    std::ostringstream s;

    s << prefix << "Read time: " << timers[READ].get_time_string() << std::endl;
    s << prefix << "Write time: " << timers[WRITE].get_time_string() << std::endl;
    s << prefix << "Database preprocess time: " << timers[PREPROCESS].get_time_string() << std::endl;
    s << prefix << "Query preprocess time: " << timers[QUERY_PREPROCESS].get_time_string() << std::endl;
    s << prefix << "Search time: " << timers[SEARCH].get_time_string() << std::endl;
    #ifdef EXTRA_TIMES
    if(timers[CANDIDATE_GENERATION].active()){
      s << prefix << "Candidate Generation time: " << timers[CANDIDATE_GENERATION].get_time_string() << std::endl;
    }
    if(timers[CANDIDATE_VERIFICATION].active()){
      s << prefix << "Candidate Verification time: " << timers[CANDIDATE_VERIFICATION].get_time_string() << std::endl;
    }
    #endif
    s << prefix << "Execution time: " << timers[EXEC].get_time_string() << std::endl;
    return s.str();
  }

  /**
   * Add up the timers from source to target.
   * @param source Source timers array.
   * @param target Target timers array.
   */
  void add_timers(
    const Timer* const source,
    Timer* const target
  ){
    for(uint8_t i=0; i < NTIMERS; ++i){
      target[i] += source[i];
    }
  }

}  // namespace timer
