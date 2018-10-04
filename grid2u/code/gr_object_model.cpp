//@header : gr_object_model.cpp
//@author : zhong

#include"include/gr_object_model.h"
#include"include/gr_constraint_generator.h"

namespace grid
{

	void ObjectModel::CollisionProcess()
	{

	}

	//��ʵ���ǵ������еķ�����
	//�Ҹо�д���������ʵ��û���õ� ��Ϊ����̬����
	//��Ϊ��������Ӧ�ò����޸��������������
	void ObjectModel::ConstraintProcess()
	{
		for (auto cg = local_cgt_.begin(); cg != local_cgt_.end(); ++cg)
		{
			(*(*cg))(*this);
		}
	}

} //namespace grid