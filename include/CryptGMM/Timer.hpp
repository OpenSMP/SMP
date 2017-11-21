//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTCONV_TIMER_HPP
#define CRYPTCONV_TIMER_HPP
#include <chrono>

using std::chrono::duration_cast;
typedef std::chrono::nanoseconds Time_t;
typedef std::chrono::high_resolution_clock Clock;
typedef Clock::duration Duration_t;
double time_as_second(const Duration_t &t) { return t.count() / 1.0e9; }
double time_as_millsecond(const Duration_t &t) { return t.count() / 1.0e6; }

#endif //CRYPTCONV_TIMER_HPP
