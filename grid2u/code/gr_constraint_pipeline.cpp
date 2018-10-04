//@source : gr_constraint_pipeline.cpp
//@author : zhong

#include"include/gr_constraint_pipeline.h"

namespace grid
{
	void ConstraintPipeline::Process()
	{
		for (auto object = constraint_omt_->begin(); object != constraint_omt_->end(); ++object)
		{
			if ((*object)) { //make sure shared ptr is invalid
				if ((*object)->CanConstraint()) {
					(*object)->ConstraintProcess();
					GlobalConstraintProcess(*(*object));
				}
			}
		}		
	}

	void ConstraintPipeline::GlobalConstraintProcess(ObjectModel& object)
	{
		for (auto cg = global_cgt_.begin(); cg != global_cgt_.end(); ++cg)
		{
			if (*cg) { //make sure shared_ptr is invalid
				(*(*cg))(object);
			} //Ӧ��������ʱ��ɾ�����Ԫ�� ĳ���Ӱ��Ч�� ������forward_list��ɾ����
		}
	}

} //namespace grid