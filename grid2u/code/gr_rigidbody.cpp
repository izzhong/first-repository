//@source : gr_ridigbody.cpp

#include"include/gr_rigidbody.h"

namespace grid
{
	void RigidBody::Update(Real duration)
	{
		Vector2 force(0, 0);
		Real moment = 0.0;
		for (auto& m : moment_)
		{
			force += m.force;
			moment += m.force.CrossProduct(m.fposition);
		}
		attr_.acceleration = force * attr_.inv_mass;
		attr_.velocity += attr_.acceleration * duration;
		attr_.position += attr_.velocity * duration;
		Real parallel = attr_.inv_rotational_inertia * attr_.inv_mass *
			(1.0 / (attr_.inv_mass + attr_.inv_rotational_inertia * (attr_.axis - attr_.position).LengthSquare()));
		attr_.angular_acceleration =
			moment * parallel;
		attr_.angular_velocity +=
			attr_.angular_acceleration * duration;
		RotateTransform(
			Vector2(0, 0),
			attr_.angular_velocity * duration, 1,
			&attr_.direction
		);
	}

} //namespace grid