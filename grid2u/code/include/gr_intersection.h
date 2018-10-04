//@header : gr_intersection.h
//@author : zhong

//TODO(zhong) : bool IntersectionTest(
//					const Geometry& geometry_1,
//					const Geometry& geometry_2,
//					bool REGION = true
//				);

//TODO(zhong) : ����ļ���ʵ�ֲ�ȡ�õ��Ǹ�ʵ���Ǹ��ķ���

//Q(zhong) : �ཻ������ʱ��ͼ�ο����ߵļ����ԵñȽϺ��� 
//				����ʱ��ͼ�ο���һ������Ҳ�ԵñȽϺ���
//				����Ӧ�����û���ѡ���ѡ�� ��Ϊһ��Ĭ��ѡ����ڵ�����������

//TODO(zhong) : ��������ͼԪ�Ĳ��Ա����϶�����ת��Ϊ�������
//				��������ǵ���ֱ�ߵ�������� ����ͨ�����ʵ��

//TODO(zhong) : ���е�ֱ�� �߶� ���� �ཻ���� �����Ի�Ϊ ��������
//				��������������㷨��಻����ʵ�ַ���

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