//@header : gr_physics_pipeline.h
//@author : zhong

#ifndef _GRID_PHYSICS_PIPELINE_H_
#define _GRID_PHYSICS_PIPELINE_H_

//TODO(zhong) : 直接绑定资源表应用的方式解决了一些问题
//				但是这么做有的时候会引起 引用了已经被释放的资源的风险
//				所以我必须借助其他方式帮助用户管理资源 从而减少这种风险发生的机会
//				所有的资源必须用智能指针创建 这样才不会出问题

#include"gr_precision.h"
#include"gr_collision.h"
#include"gr_constraint.h"

namespace grid
{
	constexpr Real PHYSICS_PIPELINE_UPDATE_DURATION = 0.16666666666667;

	class PhysicsPipeline
	{
	public:

		//物理管线在创建时需要绑定资源
		//因为其他的管线需要在创建的时候就指定资源
		//所以物理管线也必须在创建的时候指定资源
		//所以就不能让物理管线来管理资源
		//而是让其他的部分来管理资源
		//物理管线仅仅根据规则更新资源
		PhysicsPipeline(
			std::shared_ptr<ObjectModelTable> omt,
			Real duration = PHYSICS_PIPELINE_UPDATE_DURATION
		) :
			omt_(omt),
			collision_pipeline_(omt_),
			constraint_pipeline_(omt_)
		{	}

		//禁止复制和拷贝
		PhysicsPipeline(const PhysicsPipeline&) = delete;
		PhysicsPipeline& operator=(const PhysicsPipeline&) = delete;
		//谨慎的定义移动规则

		void SetUpdateDuration(Real duration)
		{
			duration_ = duration;
		}

		Real GetUpdateDuration() const
		{
			return duration_;
		}

		virtual void Process();

	protected:

		//保存一个这样的对象
		//表中的每一个智能指针都指向了一个分配在池中的对象
		//仅仅物理管线对其有读写能力
		//其余管线仅仅对其有只读权限
		std::shared_ptr<ObjectModelTable> omt_;

		//保存碰撞管线
		CollisionPipeline collision_pipeline_;

		//保存约束管线
		ConstraintPipeline constraint_pipeline_;

		//物理管线的更新时间
		Real duration_;

	}; // class PhysicsPipeline

} //namespace grid

#endif //_GRID_PHYSICS_PIPELINE_H_