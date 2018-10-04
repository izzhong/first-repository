//@header : gr_intersection.h
//@author : zhong

//TODO(zhong) : bool IntersectionTest(
//					const Geometry& geometry_1,
//					const Geometry& geometry_2,
//					bool REGION = true
//				);

//TODO(zhong) : 这个文件的实现采取用到那个实现那个的方案

//Q(zhong) : 相交测试有时将图形看作线的集合显得比较合理 
//				但有时将图形看作一个区域也显得比较合理
//				我们应该让用户来选择此选项 作为一个默认选项放在第三个参数上

//TODO(zhong) : 点与其他图元的测试本质上都可以转换为区域测试
//				其基础就是点与直线的区域测试 可以通过叉积实现

//TODO(zhong) : 所有的直线 线段 射线 相交测试 都可以化为 跨立测试
//				这与其他特殊的算法差距不大且实现方便

#ifndef _GRID_INTERSECTION_H_
#define _GRID_INTERSECTION_H_

#include"gr_geometry_algorithm.h"

namespace grid
{
	//PointIntersectionTest

	bool IntersectionTest(
		const Point& point_1,
		const Point& point_2,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Line& line,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Ray& ray,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Segment& segment,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Triangle& triangle,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Rectangle& rectangle,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Box& box,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Polygon& polygon,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Arc& arc,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Circle& circle,
		bool REGION = true
	);

	bool IntersectionTest(
		const Point& point,
		const Sector& sector,
		bool REGION = true
	);

	//LineIntersectionTest

	bool IntersectionTest(
		const Line& line,
		const Point& point,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line_1,
		const Line& line_2,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Ray& ray,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Segment& segment,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Triangle& triangle,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Rectangle& rectangle,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Box& box,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Polygon& polygon,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Arc& arc,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Sector& sector,
		bool REGION = true
	);

	bool IntersectionTest(
		const Line& line,
		const Circle& circle,
		bool REGION = true
	);

	//RayIntersectionTest

	bool IntersectionTest(
		const Ray& ray_1,
		const Ray& ray_2
	);

	bool IntersectionTest(
		const Segment& segment_1,
		const Segment& segment_2
	);

	bool IntersectionTest(
		const Triangle& triangle_1,
		const Triangle& triangle_2
	);

	bool IntersectionTest(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	);

	bool IntersectionTest(
		const Box& box_1,
		const Box& box_2
	);

	bool IntersectionTest(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	);

	bool IntersectionTest(
		const Arc& arc_1,
		const Arc& arc_2
	);

	bool IntersectionTest(
		const Sector& sector_1,
		const Sector& sector_2
	);

	bool IntersectionTest(
		const Circle& circle_1,
		const Circle& circle_2
	);

} //namespace grid

#endif //_GRID_INTERSECTION_H_