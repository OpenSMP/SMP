//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTGMM_TIMER_HPP
#define CRYPTGMM_TIMER_HPP
#include <chrono>

class AutoTimer {
public:
	using Time_t = std::chrono::nanoseconds;
	using Clock = std::chrono::high_resolution_clock;
	AutoTimer(double *ret) : ret_(ret) {
		stamp_ = Clock::now();
	}

	void reset() {
		stamp_ = Clock::now();
	}

	~AutoTimer() {
		if (ret_)
			*ret_ = (Clock::now() - stamp_).count() / 1.0e6;
	}

protected:
	double *ret_ = nullptr;
	Clock::time_point stamp_;
};
#endif //CRYPTGMM_TIMER_HPP
