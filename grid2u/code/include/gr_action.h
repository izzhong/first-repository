//@header : gr_action.h
//@author : zhong 

#ifndef _GRID_ACTION_H_
#define _GRID_ACTION_H_

#include<cassert>

namespace grid
{
	//������һ������
	//Ҫ��д�������Ļ� �͸�Լ����𲻴���
	class ActionGenerator
	{
	public:
		
		virtual enum ACTION : int;

		virtual void Action(ACTION action)
		{
			assert(action >= 0 && action < fnum_);
			(this->*act_table_[action])();
		}

		virtual int Fnum() const = 0;

		using ActionTable = void(ActionGenerator::*)();

	private:

		const int fnum_;

		ActionTable act_table_[];

	};

} //namespace grid

#endif //_GRID_ACTION_H_