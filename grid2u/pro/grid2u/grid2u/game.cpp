//@source : game.cpp

#include"game.h"
#include<iostream>
using namespace grid;
using namespace std;

namespace game
{
	//�⿴����������Լ�� �����Ƕ���
	//������Ӧ��ʹ�ú�����ʵ�ֶ�����
	void LittleCat::ConstraintAction()
	{
		cout << "��~" << endl;
	}

	void LittleCatAction::Miao()
	{
		cout << "��~" << endl;
	}

	void LittleCatAction::Jump()
	{
		cat_.AddMoment({Force(0, -10),cat_.GetAttribute().position});
	}

	void LittleCatAction::Move()
	{
		auto& attr = cat_.GetAttribute();
		cat_.AddMoment(
			{ attr.direction * 2,attr.position }
		);
	}

	void LittleCatAction::Turn()
	{
		auto& attr = cat_.GetAttribute();
		attr.direction.Invert();
	}

	void LittleCatAction::OpenLight()
	{
		(light_->action_)->Action(LightAction::OPEN_LIGHT);
	}

	void LittleCatAction::CloseLight()
	{
		(light_->action_)->Action(LightAction::CLOSE_LIGHT);
	}



	void fun()
	{
		
	}


} //namespace game