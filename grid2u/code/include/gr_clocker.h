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

	//ʱ�ӹ����߻���
	//�ṩһЩ�ӿں�״̬�趨
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
		�������ƣ�getTotalTime
		�������ܣ�����ԴӼ�ʱ��start�����ڵ���ʱ�䣨��������ͣʱ�䣩
		������������
		�������أ�float���͵�ʱ�䣬����Ϊ��λ
		*/
		Real total_time();

	protected:

		count_t start_count_ = 0;

		count_t previous_pause_count = 0;

		count_t accumulate_pause_count = 0;

	};


	//Ƶ�ʿ�����
	//ÿ������ʱ�䷵����
	//��ͣ����ֹͣ����
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

		//����tick֮��Ӧ������м����ͣʱ�� ����Ψһ������
		bool tick();

	protected:

		count_t previous_tick_count_ = 0;

		count_t current_count_ = 0;

		Real tick_delta_time_ = 0.0;

	};
}

#endif // _GRID_CLOCKER_H_