//@header : gr_constraint_generator.h
//@author : zhong

#ifndef _GRID_CONSTRAINT_GENERATOR_H_
#define _GRID_CONSTRAINT_GENERATOR_H_

//TODO(zhong) : ������������������û��Զ������Ͷ����Լ��

namespace grid
{
	class ObjectModel;

	class ConstraintGenerator
	{
	public:
		virtual void operator()(ObjectModel& object) = 0;
	};

	template<typename ObjectModel>
	class ConstraintGeneratorTemplate
	{
	public:
		virtual void operator()(ObjectModel& object) = 0;
	};

} //namespace grid

#endif //_GRID_CONSTRAINT_GENERATOR_H_