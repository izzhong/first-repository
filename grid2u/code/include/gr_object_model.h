//@header : gr_object.h
//@author : zhong

#ifndef _GRID_OBJECT_MODEL_H_
#define _GRID_OBJECT_MODEL_H_

#include"gr_constraint_generator.h"
#include"gr_precision.h"
#include<vector>
#include<forward_list>
#include<memory>

namespace grid
{
	class ObjectModel;

	//保存底层对象的智能指针
	using ObjectModelTable =
		typename std::forward_list<std::shared_ptr<ObjectModel>>; 

	//约束发生器表
	using ConstraintGeneratorTable =
		typename std::forward_list<std::shared_ptr<ConstraintGenerator>>;

	//所有对象模型需要对物理管线的支持函数都不需要由派生类进行重写
	//实际使用时需要继承这个对象 物理管线才能对其进行支持
	//按照人们的额使用习惯
	//习惯于在我们定义这个物体的时候同时定义其行为
	//而其余的约束发生器帮助其完成固有行为
	//所以约束发生器的处理与约束行为的处理应该分开
	//@date : 2018/9/22
	//对象可以具有很多行为 显然一个碰撞行为和约束行为是不够的
	//虽然我们可以为其定义大量的约束发生器 
	//但是那样显得不够灵活 所以或许应该为它实现一个动作处理器
	//用来为其定义动作和动作的发生条件 既 映射
	//这好像和约束器一样啊
	class ObjectModel
	{
	public:

		ObjectModel() {}

		//对应管线的三个状态
		bool CanConstraint() const {
			return can_constraint_;
		}

		bool CanCollision() const{
			return can_collision_;
		}

		bool CanUpdate() const {
			return can_update_;
		}

		void SetConstraintState(bool can) {
			can_constraint_ = can;
		}

		void SetCollisionState(bool can) {
			can_collision_ = can;
		}

		void SetUpdateState(bool can) {
			can_update_ = can;
		}

		void BindLocalConstraintGenerator(
			std::shared_ptr<ConstraintGenerator>& cg
		) { //倒序
			local_cgt_.push_front(cg);
		}

		void ClearLocalConstraintGeneratorTable()
		{
			local_cgt_.clear();
		}

		//碰撞处理应该由管线自己实现
		virtual void CollisionProcess();

		//约束处理也应该由管线来实现
		virtual void ConstraintProcess();

		virtual void Update(Real duration) {}

		//用户自定义的约束行为
		virtual void ConstraintAction(){}

		//发生碰撞时自定义的行为 由用户自定义
		virtual void CollisionAction() {}

	protected:

		bool can_collision_ = true;
		bool can_constraint_ = true;
		bool can_update_ = true;
		ConstraintGeneratorTable local_cgt_;//局部约束发生器
	};

} //namespace grid

#endif //_GRID_OBJECT_MODEL_H_