//@header : gr_physics_pipeline.h
//@author : zhong

#ifndef _GRID_PHYSICS_PIPELINE_H_
#define _GRID_PHYSICS_PIPELINE_H_

//TODO(zhong) : ֱ�Ӱ���Դ��Ӧ�õķ�ʽ�����һЩ����
//				������ô���е�ʱ������� �������Ѿ����ͷŵ���Դ�ķ���
//				�����ұ������������ʽ�����û�������Դ �Ӷ��������ַ��շ����Ļ���
//				���е���Դ����������ָ�봴�� �����Ų��������

#include"gr_precision.h"
#include"gr_collision.h"
#include"gr_constraint.h"

namespace grid
{
	constexpr Real PHYSICS_PIPELINE_UPDATE_DURATION = 0.16666666666667;

	class PhysicsPipeline
	{
	public:

		//��������ڴ���ʱ��Ҫ����Դ
		//��Ϊ�����Ĺ�����Ҫ�ڴ�����ʱ���ָ����Դ
		//�����������Ҳ�����ڴ�����ʱ��ָ����Դ
		//���ԾͲ��������������������Դ
		//�����������Ĳ�����������Դ
		//������߽������ݹ��������Դ
		PhysicsPipeline(
			std::shared_ptr<ObjectModelTable> omt,
			Real duration = PHYSICS_PIPELINE_UPDATE_DURATION
		) :
			omt_(omt),
			collision_pipeline_(omt_),
			constraint_pipeline_(omt_)
		{	}

		//��ֹ���ƺͿ���
		PhysicsPipeline(const PhysicsPipeline&) = delete;
		PhysicsPipeline& operator=(const PhysicsPipeline&) = delete;
		//�����Ķ����ƶ�����

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

		//����һ�������Ķ���
		//���е�ÿһ������ָ�붼ָ����һ�������ڳ��еĶ���
		//����������߶����ж�д����
		//������߽���������ֻ��Ȩ��
		std::shared_ptr<ObjectModelTable> omt_;

		//������ײ����
		CollisionPipeline collision_pipeline_;

		//����Լ������
		ConstraintPipeline constraint_pipeline_;

		//������ߵĸ���ʱ��
		Real duration_;

	}; // class PhysicsPipeline

} //namespace grid

#endif //_GRID_PHYSICS_PIPELINE_H_