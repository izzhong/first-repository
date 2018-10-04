//@header : gr_object_model.cpp
//@author : zhong

#include"include/gr_object_model.h"
#include"include/gr_constraint_generator.h"

namespace grid
{

	void ObjectModel::CollisionProcess()
	{

	}

	//其实就是调用所有的发生器
	//我感觉写成这个样其实是没有用的 因为将动态对象
	//作为函数参数应该不会修改派生对象的数据
	void ObjectModel::ConstraintProcess()
	{
		for (auto cg = local_cgt_.begin(); cg != local_cgt_.end(); ++cg)
		{
			(*(*cg))(*this);
		}
	}

} //namespace grid