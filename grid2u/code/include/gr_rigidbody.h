//gr_ridi////gbody.h
//TODO(zhong) : �Ҿ��ø�����Լ̳��ʵ�

#ifndef _GRID_RIGIDBODY_H_
#define _GRID_RIGIDBODY_H_

#include"gr_object_model.h"
#include"gr_matrix.h"

namespace grid
{
	struct RigidBodyAttribute 
	{
		Real		inv_mass;
		Vector2		position;
		Vector2		velocity;
		Vector2		acceleration;
		Vector2		axis;	//��λ��������ѡ�� ���λ�û���ȫ��λ��
		Real		inv_rotational_inertia;
		Vector2		direction;
		Real		angular_velocity;
		Real		angular_acceleration;
	};

	struct Moment
	{
		Vector2 force;
		Vector2 fposition;
	};

	class RigidBody :
		public ObjectModel
	{
	public:

		RigidBody() {}

		using Attribute = typename RigidBodyAttribute;

		RigidBody(const Attribute& attr) :
			attr_(attr), moment_()
		{	}

		virtual void Update(Real duration) override;

		//���������һ������ٶ�
		//positionΪ��������ĵ�λ��
		Vector2 LinearVelocity(
			const Vector2& position,
			const RigidBodyAttribute& attr
		) const {
			//v = w x r + vb
			return Vector2(
				-attr.angular_velocity * position.y,
				attr.angular_velocity  * position.x
			) + attr.velocity;
		}

		void AddMoment(const Moment& moment)
		{
			moment_.push_front(moment);
		}

		void ClearMoment()
		{
			moment_.clear();
		}

		Attribute& GetAttribute()
		{
			return attr_;
		}

		Attribute GetAttribute() const
		{
			return attr_;
		}

	protected:

		Attribute attr_;
		std::forward_list<Moment> moment_;

	}; //class RigidBody

} //namespace grid

#endif //_GRID_RIGIDBODY_H_