//file anno

#include"include/gr_clocker.h"
#include<cassert>
#include<Windows.h>

namespace grid
{
	ClockerBase::ClockerBase() :
		second_per_count_(0.0), state_(ClockerState::UNSTART)
	{
		count_t performance_frequency;
		QueryPerformanceFrequency(
			(LARGE_INTEGER*)&performance_frequency);
		second_per_count_ = 1.0 / Real(performance_frequency);
	}

	void ClockerBase::restart()
	{
		reset();
		start();
	}

	Timer::Timer() :ClockerBase() {  }

	void Timer::start()
	{
		if (state_ == ClockerState::UNSTART)
		{
			QueryPerformanceCounter((LARGE_INTEGER*)&start_count_);
		}
		else if (state_ == ClockerState::PAUSE)
		{
			count_t _T;
			QueryPerformanceCounter((LARGE_INTEGER*)&_T);
			accumulate_pause_count += (_T - previous_pause_count);
		}
		state_ = ClockerState::WORK;
	}

	void Timer::pause()
	{
		if (state_ == ClockerState::WORK)
		{
			QueryPerformanceCounter((LARGE_INTEGER*)&previous_pause_count);
			//_stopTime = currTime;
			state_ = ClockerState::PAUSE;
		}
	}

	void Timer::reset()
	{
		accumulate_pause_count = 0;
		state_ = ClockerState::UNSTART;
	}

	void Timer::restart()
	{
		reset();
		start();
	}

	Real Timer::total_time()
	{
		if (state_ == ClockerState::PAUSE)
			return second_per_count_ * Real(previous_pause_count - start_count_ - accumulate_pause_count);
		else if (state_ == ClockerState::WORK)
		{
			count_t _T;
			QueryPerformanceCounter((LARGE_INTEGER*)&_T);
			return second_per_count_ * Real(_T - start_count_ - accumulate_pause_count);
		}
		else
			return 0;
	}

	Ticker::Ticker(Real tick) : 
		ClockerBase(),
		tick_delta_time_(tick)
	{  
		assert(tick_delta_time_ >= 0.0);
	}

	void Ticker::set_tick(Real tick)
	{
		assert(tick >= 0);
		tick_delta_time_ = tick;
	}

	void Ticker::start()
	{
		state_ = ClockerState::WORK;
		QueryPerformanceCounter(
			(LARGE_INTEGER*)&previous_tick_count_);
	}

	void Ticker::pause()
	{
		state_ = ClockerState::PAUSE;
	}

	void Ticker::reset()
	{
		state_ = ClockerState::UNSTART;
	}

	bool Ticker::tick()
	{
		if (state_ == ClockerState::WORK)
		{
			QueryPerformanceCounter(
				(LARGE_INTEGER*)&current_count_);
			if (tick_delta_time_ <= second_per_count_ * 
				Real(current_count_ - previous_tick_count_))
			{
				previous_tick_count_ = current_count_;
				return true;
			}
			else
			{
				return false;
			}			
		}
		else
		{
			return false;
		}		
	}

}//namespace grid