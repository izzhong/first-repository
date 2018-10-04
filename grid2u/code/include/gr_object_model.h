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

	//����ײ���������ָ��
	using ObjectModelTable =
		typename std::forward_list<std::shared_ptr<ObjectModel>>; 

	//Լ����������
	using ConstraintGeneratorTable =
		typename std::forward_list<std::shared_ptr<ConstraintGenerator>>;

	//���ж���ģ����Ҫ��������ߵ�֧�ֺ���������Ҫ�������������д
	//ʵ��ʹ��ʱ��Ҫ�̳�������� ������߲��ܶ������֧��
	//�������ǵĶ�ʹ��ϰ��
	//ϰ���������Ƕ�����������ʱ��ͬʱ��������Ϊ
	//�������Լ����������������ɹ�����Ϊ
	//����Լ���������Ĵ�����Լ����Ϊ�Ĵ���Ӧ�÷ֿ�
	//@date : 2018/9/22
	//������Ծ��кܶ���Ϊ ��Ȼһ����ײ��Ϊ��Լ����Ϊ�ǲ�����
	//��Ȼ���ǿ���Ϊ�䶨�������Լ�������� 
	//���������Եò������ ���Ի���Ӧ��Ϊ��ʵ��һ������������
	//����Ϊ�䶨�嶯���Ͷ����ķ������� �� ӳ��
	//������Լ����һ����
	class ObjectModel
	{
	public:

		ObjectModel() {}

		//��Ӧ���ߵ�����״̬
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
		) { //����
			local_cgt_.push_front(cg);
		}

		void ClearLocalConstraintGeneratorTable()
		{
			local_cgt_.clear();
		}

		//��ײ����Ӧ���ɹ����Լ�ʵ��
		virtual void CollisionProcess();

		//Լ������ҲӦ���ɹ�����ʵ��
		virtual void ConstraintProcess();

		virtual void Update(Real duration) {}

		//�û��Զ����Լ����Ϊ
		virtual void ConstraintAction(){}

		//������ײʱ�Զ������Ϊ ���û��Զ���
		virtual void CollisionAction() {}

	protected:

		bool can_collision_ = true;
		bool can_constraint_ = true;
		bool can_update_ = true;
		ConstraintGeneratorTable local_cgt_;//�ֲ�Լ��������
	};

} //namespace grid

#endif //_GRID_OBJECT_MODEL_H_