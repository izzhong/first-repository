//@header : game.h
//@author : zhong
//@abstract : 整个引擎流程的初步测试

#ifndef _GRID_GAME_H_
#define _GRID_GAME_H_

#include<grid2u.h>

/*
	一个房间里有一只小猫
	还有几个开关和几个灯
	小猫可以走动 
	走到某一个开关前可以选择关掉这个灯
	或者打开这个灯
*/

namespace game
{
	//现在虽然定义了小猫咪的动作
	//但是动作需要触发 这应该怎么编写
	//小猫咪的行动有一些约束 我们应该为其编写约束器
	//碰撞管线还没有实现 所以暂时放弃
	//更新管线由物理管线完成即可

	class Light;
	//修改其约束方式有两种方法
	//直接重载其约束函数
	//或者为其编写约束发生器并绑定
	//我向给我的小猫写一个约束器 我要怎么做呢
	class LittleCat :
		public grid::RigidBody
	{
		friend class LittleCatAction;

	public:

		LittleCat(int miao_count = 0):
			miao_count_(miao_count)
		{	}

		//用户自定义的约束行为
		virtual void ConstraintAction() override;

		//发生碰撞时自定义的行为 由用户自定义
		virtual void CollisionAction() override
		{

		}
	
	protected:

		int miao_count_;

		//虽然我们绑定了动作发生器 
		//但是动作的发生需要条件
		//这些条件有的时候非常复杂
		//所以最好与动作发生器分开编写
		//动作在什么时候发生由另一部分来定义
		LittleCatAction* action_;
	};

	//这个动作代表了管线处理过程中可能发生的三种交互
	//1.无交互
	//2.自交互
	//3.与其他对象交互
	//接下来的问题是动作发生器如何与对象进行交互
	//绑定某一对象 还是绑定到某一对象
	//或者以对象为参数 不绑定到任何的对象
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

		//用户自定义的约束行为
		virtual void ConstraintAction() override
		{

		}

		//发生碰撞时自定义的行为 由用户自定义
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


	//首先我们需要给小猫咪写两个约束器
	//一个是重力约束器
	//一个是位置约束器

	//对于约束器 我们可以选择在对象构造的时候进行绑定
	//或者如果创建了大量的对象
	//可以仅仅创建一个备份 并把智能指针绑定到所有需要
	//绑定de 对象上
	class LittleCatGravityConstraint
	{
	public:

		LittleCatGravityConstraint(LittleCat& cat) :
			cat_(cat)
		{	}

		//碰到了与动作发生器一样的情况
		//约束器需要与物体进行绑定
		
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


	//但是这样约束是没有办法添加到管线中的‘
	//如果我把约束看作是一个模板可调用函数
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
