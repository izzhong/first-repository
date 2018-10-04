//@header : gr_bsp.h
//@author : zhong
//@abstract : 采用叶节点存储的BSP树

#ifndef _GRID_BSP_H_
#define _GRID_BSP_H_

#include"gr_geometry.h"
#include<cassert>
#include<initializer_list>

namespace grid
{

	//空间划分树在进行划分时
	//对象可能被划分到某一节点区域
	//也可能与某一条划分线发生碰撞

	struct SpaceNode
	{
		Line partition; //空间由直线划分
		SpaceNode* parent;
		SpaceNode* forward;
		SpaceNode* backward;
	};

	//二叉划分树就算划分完了我也不知道它属于哪个空间啊
	//所以在每个节点上
	struct Space
	{

	};

	//@class : BSP
	class BinarySpacePartition
	{
	public:

		//二叉空间树的划分由一系列直线实现
		//前向被划分到前向节点 后向被划分到后向节点
		BinarySpacePartition(std::initializer_list<Line> il)
		{
			for (auto it = il.begin(); it != il.end(); ++it)
			{
				head_ = Alloc_(*it);
			}
		}

		//我们可以按照遍历顺序创建这个树
		//但是这对用户来说不久显得有些困难吗
		//用户还得知道前序遍历的使用方法

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