//@header : gr_constraint_pipeline.h
//@author : zhong

#ifndef _GRID_CONSTRAINT_PIPELINE_H_
#define _GRID_CONSTRAINT_PIPELINE_H_

#include"gr_object_model.h"

namespace grid
{
	class ConstraintPipeline
	{
	public:

		ConstraintPipeline(std::shared_ptr<ObjectModelTable> omt) :
			constraint_omt_(omt),
			global_cgt_()
		{	}

		ConstraintPipeline(const ConstraintPipeline&) = delete;
		ConstraintPipeline& operator=(const ConstraintPipeline&) = delete;

		void BindGlobalConstraintGenerator(
			std::shared_ptr<ConstraintGenerator> global_cg
		) {
			global_cgt_.push_front(global_cg);
		}

		void ClearGlobalConstraintGenerator()
		{
			global_cgt_.clear();
		}

		virtual void Process();
		
	protected:

		virtual void GlobalConstraintProcess(ObjectModel& object);
		
	protected:

		std::shared_ptr<ObjectModelTable> constraint_omt_;
		ConstraintGeneratorTable global_cgt_;

	}; //class ConstraintPipeline

} //namespace grid

#endif //_GRID_CONSTRAINT_PIPELINE_H_