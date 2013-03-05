//============================================================================
// Name        : Timer.hpp
// Author      : Jonathan Webster & Max Revitt
// Description : 2nd Assessed Coursework for ACA 2013
//============================================================================

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <tbb/tbb.h>

// Time units of code easily with this class
class Timer
{
	private: 
		tbb::tick_count start_time, end_time;

	public:
		Timer(tbb::tick_count start_time) : start_time(start_time) {}

		void Start(tbb::tick_count now) { start_time = now; }

		void Stop(tbb::tick_count now) { end_time = now; }

		double getTime() { return (end_time-start_time).seconds(); }

		double StopWithResult(tbb::tick_count now) 
		{ 
			end_time = now;
			return (end_time-start_time).seconds();
		}
};


#endif