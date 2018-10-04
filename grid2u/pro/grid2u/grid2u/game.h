//@header : game.h
//@author : zhong
//@abstract : �����������̵ĳ�������

#ifndef _GRID_GAME_H_
#define _GRID_GAME_H_

#include<grid2u.h>

/*
	һ����������һֻСè
	���м������غͼ�����
	Сè�����߶� 
	�ߵ�ĳһ������ǰ����ѡ��ص������
	���ߴ������
*/

namespace game
{
	//������Ȼ������Сè��Ķ���
	//���Ƕ�����Ҫ���� ��Ӧ����ô��д
	//Сè����ж���һЩԼ�� ����Ӧ��Ϊ���дԼ����
	//��ײ���߻�û��ʵ�� ������ʱ����
	//���¹��������������ɼ���

	class Light;
	//�޸���Լ����ʽ�����ַ���
	//ֱ��������Լ������
	//����Ϊ���дԼ������������
	//������ҵ�Сèдһ��Լ���� ��Ҫ��ô����
	class LittleCat :
		public grid::RigidBody
	{
		friend class LittleCatAction;

	public:

		LittleCat(int miao_count = 0):
			miao_count_(miao_count)
		{	}

		//�û��Զ����Լ����Ϊ
		virtual void ConstraintAction() override;

		//������ײʱ�Զ������Ϊ ���û��Զ���
		virtual void CollisionAction() override
		{

		}
	
	protected:

		int miao_count_;

		//��Ȼ���ǰ��˶��������� 
		//���Ƕ����ķ�����Ҫ����
		//��Щ�����е�ʱ��ǳ�����
		//��������붯���������ֿ���д
		//������ʲôʱ��������һ����������
		LittleCatAction* action_;
	};

	//������������˹��ߴ�������п��ܷ��������ֽ���
	//1.�޽���
	//2.�Խ���
	//3.���������󽻻�
	//�������������Ƕ�������������������н���
	//��ĳһ���� ���ǰ󶨵�ĳһ����
	//�����Զ���Ϊ���� ���󶨵��κεĶ���
	class LittleCatAction
	{
	public:

		LittleCatAction(LittleCat& cat) :
			cat_(cat)
		{	}

		enum ACTION : int
		{
			MIAO,
			JUMP,
			MOVE,
			TURN,
			OPEN_LIGHT,
			CLOSE_LIGHT
		};

		using FuntionTable = void(LittleCatAction::*)();

		int Fnum() const
		{
			return fnum_;
		}

		void Action(ACTION action)
		{
			assert(action >= 0 && action < fnum_);
			(this->*ftable_[action])();
		}

	private:

		void Miao();

		void Jump();

		void Move();

		void Turn();

		void OpenLight();

		void CloseLight();

	private:

		LittleCat & cat_;

		Light* light_;

		int fnum_ = 6;

		FuntionTable ftable_[6] = {
			&LittleCatAction::Miao,
			&LittleCatAction::Jump,
			&LittleCatAction::Move,
			&LittleCatAction::Turn,
			&LittleCatAction::OpenLight,
			&LittleCatAction::CloseLight
		};
	};

	class Light :
		public grid::Particle
	{
		friend class LightAction;

	public:

		//�û��Զ����Լ����Ϊ
		virtual void ConstraintAction() override
		{

		}

		//������ײʱ�Զ������Ϊ ���û��Զ���
		virtual void CollisionAction() override
		{

		}

		LightAction* action_;

	private:

		bool state;		
	};

	class LightAction
	{
	public:

		LightAction(Light& light) :
			light_(light)
		{	}

		enum ACTION : int
		{
			OPEN_LIGHT,
			CLOSE_LIGHT
		};



		using FunctionTable = void(LightAction::*)();

		int Fnum() const
		{
			return fnum_;
		}

		void Action(ACTION action)
		{
			assert(action >= 0 && action < fnum_);
			(this->*act_table_[action])();
		}

	private:

		void OpenLight()
		{
			light_.state = true;
		}

		void CloseLight()
		{
			light_.state = false;
		}

		Light & light_;

		int fnum_ = 2;

		FunctionTable act_table_[2] = {
			&LightAction::OpenLight,
			&LightAction::CloseLight
		};
	};


	//����������Ҫ��Сè��д����Լ����
	//һ��������Լ����
	//һ����λ��Լ����

	//����Լ���� ���ǿ���ѡ���ڶ������ʱ����а�
	//������������˴����Ķ���
	//���Խ�������һ������ ��������ָ��󶨵�������Ҫ
	//��de ������
	class LittleCatGravityConstraint
	{
	public:

		LittleCatGravityConstraint(LittleCat& cat) :
			cat_(cat)
		{	}

		//�������붯��������һ�������
		//Լ������Ҫ��������а�
		
		void operator()()
		{
			cat_.AddMoment(
				{Force(0,9.8),cat_.GetAttribute().position}
			);
		}

	private:

		LittleCat & cat_;
	};

	class LittleCatPositionConstraint
	{
	public:

		LittleCatPositionConstraint(LittleCat& cat) :
			cat_(cat)
		{	}

		void operator()()
		{
			auto& attr = cat_.GetAttribute();
			if (attr.position.y > 600) {
				attr.position.y = 600;
			}
			if (attr.position.x < 0) {
				attr.position.x = 0;
			}
			if (attr.position.x > 1067) {
				attr.position.x = 1067;
			}
		}

	private:

		LittleCat & cat_;
	};


	//��������Լ����û�а취��ӵ������еġ�
	//����Ұ�Լ��������һ��ģ��ɵ��ú���
	template<typename ObjectModel>
	class ConstraintGenerator
	{
	public:
		virtual void operator()(ObjectModel& object) = 0;
	};

	class LittleCatConstraintGenerator :
		public ConstraintGenerator<LittleCat>
	{	};

	class LittleCatGravityConstraint :
		public LittleCatConstraintGenerator
	{
	public:
		virtual void operator()(LittleCat& cat) override
		{
			cat.AddMoment(
				{ Force(0,9.8),cat.GetAttribute().position }
			);		
		}
	};

} //namespace game

#endif //_GRID_GAME_H_
