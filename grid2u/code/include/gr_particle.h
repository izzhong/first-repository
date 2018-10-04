//@header : gr_particle.h
//@author : zhong

#ifndef _GRID_PARTICLE_H_
#define _GRID_PARTICLE_H_

#include"gr_object_model.h"
#include"gr_vector.h"

namespace grid
{
	struct ParticleAttribute
	{
		Real inv_mass;
		Vector2 position;
		Vector2 velocity;
		Vector2 acceleration;
	};

	using Force = Vector2;

	//class Particle
	class Particle :
		public ObjectModel
	{
	public:

		using Attribute = typename ParticleAttribute;

		Particle() {}

		Particle(const Attribute& attr) :
			attr_(attr),force_()
		{	}

		virtual void Update(Real duration) override
		{
			attr_.acceleration = force_ * attr_.inv_mass;
			attr_.velocity += attr_.acceleration * duration;
			attr_.position += attr_.velocity * duration;
		}

		Attribute& GetAttribute()
		{
			return attr_;
		}

		void AddForce(const Force& force)
		{
			force_ += force;
		}

		void ClearForce()
		{
			force_.Clear();
		}

	private:
		
		Attribute attr_;
		Force force_;

	}; //class Particle

} //namespace gid

#endif //_GRID_PARTICLE_H_