//@source : gr_physics_pipeline.cpp
//@author : zhong

#include"include/gr_physics_pipeline.h"
#include<iostream>

namespace grid
{
	void PhysicsPipeline::Process()
	{		
		if (omt_) {
			constraint_pipeline_.Process();
			collision_pipeline_.Process();
			for (auto object = omt_->begin(); object != omt_->end(); ++object)
			{
				if ((*object)) {
					(*object)->Update(duration_);
				}
			}
		}
	}

} //namespace grid