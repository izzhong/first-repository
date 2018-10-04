//file anno

#ifndef _GRID_CLOCKER_H_
#define _GRID_CLOCKER_H_

#include"gr_precision.h"

namespace grid
{
	//
	enum class ClockerState : int
	{
		UNSTART,
		WORK,
		PAUSE
	};

	//时钟管理者基类
	//提供一些接口和状态设定
	class ClockerBase
	{
	public:

		ClockerBase();

		virtual void start() = 0;

		virtual void pause() = 0;

		virtual void reset() = 0;

		virtual void restart();

		virtual ClockerState get_state() const { return state_; }

		bool is_unstart() const { return state_ == ClockerState::UNSTART; }

		bool is_pause() const { return state_ == ClockerState::PAUSE; }

		bool is_work() const { return state_ == ClockerState::WORK; }

	protected:

		Real second_per_count_;

		ClockerState state_;

	protected:

		using count_t = unsigned long long;
	}; //class ClockerBase


	class Timer :
		public ClockerBase
	{
	public:

		Timer();

		virtual void start() override;

		virtual void pause() override;

		virtual void reset() override;

		virtual void restart() override;

		/*
		函数名称：getTotalTime
		函数功能：获得自从计时器start后到现在的总时间（不包括暂停时间）
		函数参数：无
		函数返回：float类型的时间，以秒为单位
		*/
		Real total_time();

	protected:

		count_t start_count_ = 0;

		count_t previous_pause_count = 0;

		count_t accumulate_pause_count = 0;

	};


	//频率控制器
	//每隔设置时间返回真
	//暂停可以停止工作
	class Ticker :
		public ClockerBase
	{
	public:

		Ticker(Real tick);

		void start();

		void pause();

		void reset();

		void set_tick(Real tick);

		Real get_tick() const { return tick_delta_time_; }

		//两次tick之间应该清空中间的暂停时间 这是唯一的区别
		bool tick();

	protected:

		count_t previous_tick_count_ = 0;

		count_t current_count_ = 0;

		Real tick_delta_time_ = 0.0;

	};
}

#endif // _GRID_CLOCKER_H_