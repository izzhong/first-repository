//@header : gr_collision_pipeline.h
//@author : zhong

#ifndef _GRID_COLLISION_PIPELINE_H_
#define _GRID_COLLISION_PIPELINE_H_

#include"gr_object_model.h"

namespace grid
{

	//��Ҫ������ײ��
	//��ײ������ײ�����
	struct CollisionObjectPair
	{
		std::shared_ptr<ObjectModel> first;
		std::shared_ptr<ObjectModel> second;
	};

	//��ײ������ײ�����
	using CollisionObjectTable =
		typename std::forward_list<CollisionObjectPair>;

	class CollisionPipeline
	{
	public:

		CollisionPipeline(std::shared_ptr<ObjectModelTable> omt) :
			omt_(omt)
		{	}

		//�����������Ҫ���ݶ����������ײ��
		//������ײ��֮�� 
		//�����ɵ���ײ��Ķ���ִ����ײ����
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

		//ÿ��process֮�󱣴���ײ��
		CollisionObjectTable collision_omt_;

		//��ײ������Ҫ�����ض����δ���BSP��
		//��BSP��������һ���ĵ��� 
		//��BSP����ÿһ��Ҷ�ڵ��ϱ������İ�Χ��νṹ

	};

} //namespace grid

#endif //_GRID_COLLISION_PIPELINE_H_