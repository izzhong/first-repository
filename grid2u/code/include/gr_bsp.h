//@header : gr_bsp.h
//@author : zhong
//@abstract : ����Ҷ�ڵ�洢��BSP��

#ifndef _GRID_BSP_H_
#define _GRID_BSP_H_

#include"gr_geometry.h"
#include<cassert>
#include<initializer_list>

namespace grid
{

	//�ռ仮�����ڽ��л���ʱ
	//������ܱ����ֵ�ĳһ�ڵ�����
	//Ҳ������ĳһ�������߷�����ײ

	struct SpaceNode
	{
		Line partition; //�ռ���ֱ�߻���
		SpaceNode* parent;
		SpaceNode* forward;
		SpaceNode* backward;
	};

	//���滮�������㻮��������Ҳ��֪���������ĸ��ռ䰡
	//������ÿ���ڵ���
	struct Space
	{

	};

	//@class : BSP
	class BinarySpacePartition
	{
	public:

		//����ռ����Ļ�����һϵ��ֱ��ʵ��
		//ǰ�򱻻��ֵ�ǰ��ڵ� ���򱻻��ֵ�����ڵ�
		BinarySpacePartition(std::initializer_list<Line> il)
		{
			for (auto it = il.begin(); it != il.end(); ++it)
			{
				head_ = Alloc_(*it);
			}
		}

		//���ǿ��԰��ձ���˳�򴴽������
		//��������û���˵�����Ե���Щ������
		//�û�����֪��ǰ�������ʹ�÷���

	private:

		SpaceNode * Alloc_(const Line& line)
		{
			SpaceNode* node = nullptr;
			node = new SpaceNode{ line,nullptr,nullptr,nullptr };
			assert(node);
			return node;
		}

		SpaceNode * head_;

	};

} //namespace grid

#endif //_GRID_BSP_H_