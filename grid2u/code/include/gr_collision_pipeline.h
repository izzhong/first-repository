//@header : gr_collision_pipeline.h
//@author : zhong

#ifndef _GRID_COLLISION_PIPELINE_H_
#define _GRID_COLLISION_PIPELINE_H_

#include"gr_object_model.h"

namespace grid
{

	//需要定义碰撞表
	//碰撞表由碰撞对组成
	struct CollisionObjectPair
	{
		std::shared_ptr<ObjectModel> first;
		std::shared_ptr<ObjectModel> second;
	};

	//碰撞表由碰撞对组成
	using CollisionObjectTable =
		typename std::forward_list<CollisionObjectPair>;

	class CollisionPipeline
	{
	public:

		CollisionPipeline(std::shared_ptr<ObjectModelTable> omt) :
			omt_(omt)
		{	}

		//这个函数就需要根据对象表生成碰撞表
		//生成碰撞表之后 
		//对生成的碰撞表的对象执行碰撞处理
		virtual void Process()
		{
			collision_omt_.clear();
		}

		CollisionObjectTable& GetCollisionObjectTable()
		{
			return collision_omt_;
		}

	protected:

		std::shared_ptr<ObjectModelTable> omt_;

		//每次process之后保存碰撞表
		CollisionObjectTable collision_omt_;

		//碰撞管线需要根据特定地形创建BSP树
		//此BSP树代表着一定的地形 
		//再BSP树的每一个叶节点上保存对象的包围层次结构

	};

} //namespace grid

#endif //_GRID_COLLISION_PIPELINE_H_