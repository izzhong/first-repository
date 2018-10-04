//@header : gr_grometry_algorithm.h
//@author : zhong
//@date	  : 2018/9

//TODO(zhong) : �ཻ����д����ͬ�����������ڹ��ߵ�ģ�廯
//				bool IntersectionTest(const Geometry left,const Geometry& right);

#ifndef _GRID_GEOMETRY_ALGORITHM_H_
#define _GRID_GEOMETRY_ALGORITHM_H_

#include"gr_geometry.h"

namespace grid
{
	
	//------ֱ���߶������㷨----------

	void TwoPointsLineEquation(
		const Point& A, const Point& B,
		Real* coeA, Real* coeB, Real* coeC
	);

	void PointNormalLineEquation(
		const Point& P,
		const Vector2& normal,
		Real* coeA, Real* coeB, Real* coeC
	);


	//��ʵ������߶ο�������
	bool SegmentStraddeleTest(
		const Segment& segment_1,
		const Segment& segment_2
	);

	////��A B��ֱ��CD���� ������
	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Point& LineC,
		const Point& LineD
	);

	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Line& line
	);

	//�ж����� OC �Ƿ������� OA OB ֮��
	bool VectorRegionTest(
		const Vector2& OA, 
		const Vector2& OB,
		const Vector2& OC
	);

	//-----��Զ���㷨------

	////����������������Զ����±�
	int GetFartestPointIndex(
		const Vector2 direction,
		int size, const Point* point_array
	);

	Point GetFartestPoint(
		const Vector2 direction,
		int size, const Point* point_array
	);

	PointNode* GetFartestPointPtr(
		const Vector2& direction,
		const PointSet& pointset
	);

	Point GetFartestPoint(
		const Vector2& direction,
		const PointSet& pointset
	);

	//-----�����------

	Point ClosestPointPtSeg(
		const Point& point, const Segment& segment
	);

	Point ClosestPointPtRec(
		const Point& point,const Rectangle& rectangle
	);

	Point ClosestPointPtTri(
		const Point& point, const Triangle& triangle
	);

	Point ClosestPointPtBox(
		const Point& point, const Box& box
	);

	Point ClosestPointPtPol(
		const Point& point, const Polygon& polygon
	);

	Point ClosestPointSegSeg(
		const Segment& segment_1, const Segment& segment_2
	);

	Point ClosestPointRayRay(
		const Ray& ray_1, const Ray& ray_2
	);

	Point ClosestPointSegTri(
		const Segment& segment, const Triangle& triangle
	);

	Point closestPointTriTri(
		const Triangle& triangle_1, const Triangle& triangle_2
	);

	//-----�����㷨-------

	////@abstract : calculate the distance between two points
	Real CalculatePointsDistance(
		const Point& point_1,
		const Point& point_2
	);

	//@abstract : calculate the cloest distance between 
	//			  a point and a segment
	Real CalculatePointToSegmentDistance(
		const Point& point,
		const Segment& segment
	);

	//----ͶӰ�����㷨-------
	//@abstract : ����ͶӰ����
	Real CalculateProjectionLength(
		const Vector2& axis,
		int size, const Point* point_array,
		Real* min, Real* max,
		int* min_nptr = nullptr, int* max_nptr = nullptr
	);

	//�����������Ը���ָ�� ��ָ��Ͳ�������Ϳ�����
	Real CalculateProjectionLength(
		const Vector2& axis,
		const PointSet& pointset,
		Real* min, Real *max,
		PointNode* min_nptr = nullptr, PointNode* max_nptr = nullptr
	);

	Real CalculateProjectionLength(
		const Vector2& axis,
		const Circle& circle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Arc& arc,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Sector& sector,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Triangle& triangle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Rectangle& rectangle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Box& box,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Polygon& polygon,
		Real* min, Real*max
	);

	//��Χ���㷨
	Rectangle RectangleBounding(const Polygon& polygon);

	//--------���������㷨-------

	void FindHull(
		const PointSet& Sk,
		PointNode* P, PointNode* Q,
		PointSet* polygon
	);

	void QuickHull(
		int size, Point* point_array,
		PointSet* pointset
	);

	////ʵ��͹����μ���㷨
	///*
	//�ǳ���
	//͹����ξ���������һ���߻�һ���� ʣ�µ����е㶼�������ߵ�һ��
	//�����жϵ����ߵķ����ϵ �õ�� ���� ��˶�����
	//�����һ���ǳ���Ҫ�����ǿ���ͨ�����ķ����ж���ʸ���໥֮���˳��ʱ���ϵ
	//�㷨
	//˳��ѡȡ������ �����������ȷ���������ε���ת����
	//ʣ�µ���������ѡȡ�� һ����ת�������ʼ����ͬ
	//��������һ���������
	//�������͹�����
	//�㷨ʱ�临�Ӷ� o[n]
	//*/
	bool IsConvexPolygon(
		const Polygon& polygon
	);

	//------------------------------------------------


	//�������㷨

	template<typename Geometry1,typename Geometry2>
	bool SeparationDetection(
		const Geometry1& gry_1, const Geometry& gry_2
	);

	bool CircleSeparationDetection(
		const Circle& circle_1,
		const Circle& circle_2
	);

	bool RectangleSeparationDetection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	);

	bool BoxSeparationDetection(
		const Box& box_1,
		const Box& box_2
	);

	bool PolygonSeparationDetection(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	);

	//�ཻ���

	bool IntersectionDetection();

	//------��ײ����㷨----

	//�����ཻ״̬
	enum class GeometryIntersectionState : __int8
	{
		SEPERATION,
		OUTER_CONTACT,
		INTERSECTION,
		INTER_CONTACT,
		CONTAINT
	};

	struct CollisionState
	{
		Vector2 collision_vector;
		Real collision_depth;
		Point closest_point[2];
	};

	//����ѧ�ཻ����㷨���
	//û�б�Ҫÿ����ײ�������˸���
	//�󲿷�ʱ������ֻ����ֱ��������û����ײ
	struct IntersectionResult
	{
		GeometryIntersectionState state;
		Vector2	intersection_vector;
		Real	intersection_depth;
		Point	closet_point[2];//����������ײ��		
	};

	bool IsPointSegmentCollineation(
		const Point& point,
		const Segment& segment
	);

	bool IsPointOnSegment(
		const Point& point,
		const Segment& segment
	);

	bool IsPointInTriangle(
		const Point& point,
		const Triangle& triangle
	);

	bool IsPointOnTriangle(
		const Point& point,
		const Triangle& triangle
	);

	bool IsPointOutTriangle(
		const Point& point,
		const Triangle& triangle
	);

	////rectangle
	///*
	//������ת�ľ��ε���ײ���
	bool IsRectangleIntersection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	);

	//point to point intersection detection
	void PointToPointIntersectionDetection(
		const Point& point_1,
		const Point& point_2,
		IntersectionResult *result
	);

	//�Լ򵥼���ͼ��ֱ��Ӧ��EPA�㷨 
	//�������ʽֱ�Ӽ���Ϳ����� 
	//����Щ��������ʲô�Ŀ����
	//��������߶ε���̾��� �� ��̾����
	void PointToSegmentIntersectionDetection(
		const Point& point,
		const Segment& segment,
		IntersectionResult* result
	);

	void PointToTriangleIntersectionDetection(
		const Point& point,
		const Triangle& triangle,
		IntersectionResult* result
	);

	void PointToCircleIntersectionDetection(
		const Point& point,
		const Circle& circle,
		IntersectionResult* result
	);
	
	void PointToRectangleIntersectionDetection(
		const Point& point,
		const Rectangle& rectangle,
		IntersectionResult* result
	);

	void PointToBoxIntersectionDetection(
		const Point& point,
		const Box& box,
		IntersectionResult* result
	);

	void PointToSectorIntersectionDetection(
		const Point& point,
		const Sector& sector,
		IntersectionResult* result
	);

	void PointToArcIntersectionDetection(
		const Point& point,
		const Arc& arc,
		IntersectionResult* result
	);

	////�жϵ��Ƿ��ڶ������
	//inline bool point_in_polygon(Vector2 _point, polygon _pol)
	//{
	//	//���߷�
	//	//������Զ��ȡһ�� ����ָ�����ƶ�
	//	//�����������ε�һ���������������� 
	//	//���������һ�� ��ʾ���˶����

	//	//��֮�����������������ڶ������
	//	//����ż�������ʾ���ڶ������

	//	//��ʵ����ʵ��������ž�������ָ����λһ�� ��������λ��������ֱ�����������бߵĽ����

	//	//���ĳЩ������� �ڶԶ���ε�ÿ���߽����ж�ʱ Ҫ����һ���������
	//	/*
	//	���赱ǰ����ı���p1p2
	//	1.�������p1p2���� ��ֱ���ж���p�ڶ������
	//	2.�����p�������������ô���p1��p2 ��ô�������ᱻ��������
	//	����������Ĵ���ԭ����: ���p��y������p1,p2�н�С��y������ͬ ������������
	//	3.�����p������������p1 p2 ƽ�� �����������
	//	*/

	//	/*
	//	�����Ĺ��̴��� ��������һ������ε�������ο����ж� �Բ��ٶ������Χ������������ų�
	//	�����ε���������㷨
	//	��������ε����нڵ��ҳ����������ϵ���ֵ����
	//	*/
	//}
	void PointToPolygonIntersectionDetection(
		const Point& point,
		const Polygon& polygon,
		IntersectionResult* result
	);

	////�ж��߶��Ƿ��ཻ
	//�������������
	//1.�����ų�����
	//�����߶�p1p2Ϊ�Խ��ߵľ���ΪR1,�����߶�Q1Q2Ϊ�Խ��ߵľ���ΪR2
	//���R1 R2���ཻ �����߶β����н���
	//2.��������
	//������߶��ཻ �����߶α�Ȼ�໥�����Է�
	//��ν���� ָ����һ���߶ε����ϵ�ֱ�λ����һ�߶�����ֱ�ߵ�����
	//��p1p2����q1q2 ��ʸ�� p1-q1 �� p2-q2λ��ʸ��q2-q1������
	//������
	//(p1-q1)x(q2-q1) * (p2-q1)x(q2-q1) < 0
	void SegmentIntersectionDetection(
		const Segment& segment_1,
		const Segment& segment_2,
		IntersectionResult* result
	);

	void SegmentRectangleIntersectionDetection(
		const Segment& segment,
		const Rectangle rectangle,
		IntersectionResult* result
	);

	//���������ܽ�
	//�߶�������ͼ�ε��ཻ��������������
	/*
		����ִ�п��ټ�� ���μ��
		Ȼ��ִ�п������
	*/

	////ֱ�߼��
	////���ߺ�Բ�εļ��
	//inline bool line_circle(
	//	Vector2 _line_v, Vector2 _point, circle _cir
	//)
	//{
	//	//�жϵ㵽ֱ�ߵľ�����뾶�Ĳ�ֵ

	//	//������Ҫ�ж� �Ƿ�ֱ���ཻʱ �������߲����ཻ 
	//	//ֻ��Ҫ���� �ӷ���㵽Բ�ĵĵ�������ֱ����������������0�Ĺ�ϵ����
	//}
	void SegmentCircleIntersectionDetection(
		const Segment& segment,
		const Circle& circle,
		IntersectionResult* result
	);

	
	//Բ�κ;���
	//inline bool circle_rectangle(
	//	cir* _cir, rectangle _rec
	//)const
	//{
	//	Vector2 close;
	//	//����ȷ��Բ�������������ĵ�
	//	if (_cir->x <= _rec->x) close.x = _rec->x;
	//	else if (_cir->x >= _rec->x + _rec->width) close.x = _rec->x + _rec->width;
	//	else close.x = _cir->x;
	//	//yҲһ��

	//	//Ȼ�����˵㵽Բ�ĵľ��벢��뾶�Ƚ�
	//	return (close - Vector2(_cir->x, _cir->y)).lengthSquare() < _cir->radius * _cir->radius;
	//}
	void RectangleCircleIntersectionDetection(
		const Rectangle& rectangle,
		const Circle& circle,
		IntersectionResult* result
	);

	////Բ������ת����
	///*
	//��Ȼ������������Ϊ���������ת �����ж�����Բ���Ƿ�������ײ�����ϻ����ҳ���������Բ�ĵ������
	//���ǿ��԰�����ͼ�ο����Ի���������������ת
	//ʹ֮��ɳ�����ж�
	//*/
	void BoxCircleIntersectionDetection(
		const Box& box,
		const Circle& circle,
		IntersectionResult* result
	);

	//��������ת��
	void RectangleBoxIntersectionDetection(
		const Rectangle& rectangle,
		const Box& box,
		IntersectionResult* result
	);

	//��ת����ײ���
	void BoxIntersectionDetection(
		const Box& box_1,
		const Box& box_2,
		IntersectionResult* result
	);

	void CircleIntersectionDetection(
		const Circle& circle_1,
		const Circle& circle_2,
		IntersectionResult* result
	);



	////���һ��������ײ���㷨
	////���κ�Բ�ε���ײ�ж�
	//bool sector_circle()
	//{
	//	Vector a(sector.vertex[0] - sector.center);
	//	Vector b(sector.vertex[1] - sector.center);
	//	Vector c(circle.center - sector.center);
	//	//�ж�Բ���Ƿ�������������
	//	if (a x c * c x b > 0 && a x c * a x b > 0)
	//	{
	//		//ֱ��ʹ��Բ����ײ����㷨
	//	}
	//	else
	//	{
	//		//����Բ�����߶ε���̾���
	//		//��뾶���Ƚϼ���
	//	}
	//}
	//Բ�κ�����
	//��������ļ��ʮ�ָ��� ����������������ô��
	//��ʱ����ʵ��
	void CircleSectorIntersectionDetection(
		const Circle& circle,
		const Sector& sector,
		IntersectionResult* result
	);


	////SAT
	bool SAT(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	);


	Vector2 GetClosetPointOnPolygon(
		const Point& point,
		const Polygon& polygon
	);



	////Բ�������ε���ײ���
	///*
	//ֻ��Ҫ��Բ��Ͷ�䵽һ��ͶӰ���ϼ���
	//���������Բ�������ζ����������һ�������
	//*/
	bool CirclePolygonSAT(
		const Circle& circle,
		const Polygon& polygon
	);



		////-------------------------------------------------------------------
		///*
		//�ǳ�ţ�Ƶ�һ��������
		//���ڼ����ཻ���
		//GJK�㷨

		//!!:���������״�ص������ཻ��ô��������״���ɿƷ�˹��������ĵ����Χ��״�����ԭ��
		//*/

		////GJK�㷨֧�ֺ���
		//Vector2 support(polygon _poy1, polygon _poy2, Vector2 _dir)
		//{
		//	//�������������ڸ�����������Զ�ĵ�
		//	//��ʵ�����ҵ��������������������������Զ��������
		//	//��Ϊһ�����Ʒ�˹����
		//	//�� ���������Ϳ�����

		//	//��ô��ô���������������������ϵ���Զ�ĵ���??
		//	for (size_t i = 0; i < _poy1.side; ++i)
		//	{
		//		//������������ڸ��������ϵķ��� ȡ�����ֵ
		//		//�ڼ�������� �������ֵ�ڵ���±�
		//	}
		//	//����poy2 ��ʹ��-dir���м��� ȡ��������
		//	for (size_t i = 0; i < _poy2.side; ++i)
		//	{
		//		//������������ڸ��������ϵķ��� ȡ�����ֵ
		//		//�ڼ�������� �������ֵ�ڵ���±�
		//	}

		//	//����������������Ʒ�˹����

		//}

		////��������
		//void GJK()
		//{
		//	vector d = ; //����ķ���ѡ��Ҳ�Ե�����������һЩӰ��
		//				 //������ε����Ŀ����ʱ ���ǿ��Բ�ȡһЩ����
		//				 //ʲô�����ܴ���һ��������� ��Ȼ������!!����ʾ��������ʹ�������������������Ϊ����յ��������Ϊ��ʼ����
		//				 //���ܾ���Ӧ���Ǵ�ֱ����
		//				 ///�������������غ� ���ľ�һ�� ����ֱ���ж��ཻ
		//				 //��Ҫ��һ��
		//	simplex.add(support(A, B, d));
		//	d.negate();
		//	while (true)
		//	{
		//		simplex.add(support(A, B, d));
		//		//���ո���ӵ��������Ƿ�����ԭ�� 
		//		if (simplex.getLast().dot(d) <= 0)
		//		{
		//			//�ж� �������� �� d ���ɵ����Ʒ�˹�����
		//			//�Ƿ�����������Ͽ����ԭ��
		//			//���Ǵ�ԭ�� ��d��ֱ�ķ���һ��ֱ��
		//			//��������ԭ�� ��ôָ�򽫰������������߷ֿ� 
		//			//Ҳ���� ��� >0
		//			//������Ϊ0 ��ô�����ɵĵ����ֱ����
		//			//������<0 ��ô�����ɵĵ��û�������������ԭ�� 
		//			//���ﻹ��Ҫ�о� == 0������ʲô��� �ǽӴ���??


		//			return false;
		//		}

		//		else
		//		{
		//			//����������Ҫȷ��ԭ���Ƿ��ڵ�������
		//			if (containsOrigin(simplex, d))
		//			{
		//				//�����������ײ
		//				return true;
		//			}

		//		}
		//	}
		//}

		////�ڼ���ʱ  
		////�������ε�Ӵ� ��ô����һ���������ڼ���ʱ������������
		////
		////���кܶ������Ҫ������ �����ı߽Ӵ�  ����  ������Ӵ�  �����߽Ӵ� ��������

		////�����Ż� 
		////������Ǽ������һ�������εĵ� ����ԭ�� ��ô���ǿ϶���ײ�� �Ӵ���ײ
		////�������������߰���ԭ�� ��ôҲ�϶���ײ�� 


		//bool containsOrigin(simplex s, Vector d)
		//{
		//	//һ�㶼��ȡ�ո���ӵĵ�Ϊ���

		//	{
		//		//���㷽���ʱ����Ҫ����ԭ���Զ�ĵ�ȥ��
		//		//ÿ�μ��������Ҫ������

		//		//���ǻ�Ҫ�������ӵ�����
		//		//d.getDirection(simplex);
		//		//ѡ����һ����������
		//		//Vector2 d;
		//		//Vector2 a = support(.., .., d);
		//		//Vector2 b = support(.., .., -d);
		//		//Vector2 AB = a - b;
		//		//Vector2 A0 = a - Vector2(0, 0);
		//		//d = AB X AO X AB; //��һ�仰�ҵ��Լ������ ����Ȼ�Ƕ�ά�� �������յĽ��Ҳȷʵ���ŵ�������ƽ��
		//		// a x b x c = b(a*c) - a(b*c)
		//		//��ʵ����ֻ�����˸ո�����������������ɵ�������һ����ֱ����
		//		//�Ѳ��ɷ����н��� �������������õ�

		//		//!!����
		//		//�����������ʵ��ͨ��������ѡȡ����һ��������
		//		//A B ��˳��������ν��
		//		//����һ��Ҫ���������ļ�����ʽ
		//		//BA X BO X BA
		//		//Ҳ����
		//		//��������һ��ѡ�������ָ��ԭ����һ��

		//	}


		//	//һ��Ҫ��øո���ӵĵ� ��Ϊa
		//	//��������Ϊb c ˳���޹�
		//	//����Ҫ�жϵ��� ԭ�㵽���� ab ���� ������ ac���� ����������������
		//	//ͨ�� acxabxab���Ի�� ab ָ������ķ�����
		//	//ͨ�� abxacxac���Ի�� ac ָ������ķ�����
		//	//���Ƕ��п��ܵ���0�� ���ԭ��λ�ڵ������ϵ�����ζ��ʲô???
		//	//��������ӵĵ�
		//	a = simplex.getLast();
		//	//compute AO
		//	ao = a.negate();
		//	if (simplex.points.size() == 3)
		//	{
		//		//then its the triangle case
		//		//get b and c
		//		b = simplex.getB();
		//		c = simplex.getC();
		//		//compute the edges
		//		ab = b - c;
		//		ac = c - a;
		//		//compute the normals
		//		abPrep = tripleProduct(ac, ab, ab);
		//		acPrep = tripleProduct(ab, ac, ac);
		//		//is the origin in R4
		//		if (abPrep.dot(ao) > 0)
		//		{
		//			//remove point c
		//			simplex.remove(c);
		//			//set the new direction to abPrep
		//			d.set(abPrep);
		//		}
		//		else
		//		{
		//			//is origin in R3
		//			if (acPrep.dot(ao) > 0)
		//			{
		//				//remove point b
		//				simplex.remove(b);
		//				//set the new direction to acPrep
		//				d.set(acPrep);
		//			}
		//			else
		//			{
		//				//otherwise we hnow its in R5 wo we can return 
		//				return true;
		//			}
		//		}
		//	}
		//	else
		//	{
		//		//then it is the line segment case
		//		b = simplex.getB();
		//		//compute AB
		//		ab = b - a;
		//		//get the perp to AB in the direction of the origin
		//		abPerp = tripleProduct(ab, ab, ab);

		//		//��������һ�����
		//		//��ͨ��������ѡ������һ������ ���������ָ��ԭ��ķ���ʱ
		//		//�������������������
		//		//��ʱ��ζ�� �߶θպ�ͨ��ԭ��
		//		/*
		//		The catch here is what happens when O lies on the line?
		//		If that happens the perp will be a zero vector and will cause the check on line 11 to fail.
		//		This can happen in two places: 1) inside the Minkowski Sum and 2) on the edge of the Minkowski Sum.
		//		The latter case indicates a touching contact rather than penetration so you will need to make a decision on whether this is considered a collision or not.
		//		In either case, you can use either the left or right hand normal of AB as the new direction.
		//		*/
		//		//�п��ܷ����Ӵ���ײ


		//		//set the direction to abPerp
		//		d.set(abPerp);
		//	}
		//}


		////---------------------------------�������������-----------------------------------------

		////��������֮�����̾���������ǵ����ɷ�˹������״��ԭ����������
		///*
		//�������͹������ɷ�˹������״��û�а���ԭ�� �����������岻���ཻ
		//��� ���ǲ��ò�ȡ�����ķ������ϼ����Χԭ��ĵ�����
		//�����뷨�õ�һ����ԭ������ĵ�����
		//����ĵ��������������ɷ�˹����ı߽���
		//*/

		///*
		//��һ���ֵĴ��뿴������GJK����
		//��ͬ���� ���ۺ������Ƕ�ֻ����2���� ��ά�������������
		//�����ڵ�������Ѱ����ԭ������ĵ� ������Ѱ��ԭ������Voronoi����

		//�ⲿ�ִ����������͹�����
		//���� ���������ཻ ��ô�㷨Ҳ�����
		//����һ��û������ ��Ϊ����һ�㶼���ȼ����ײ

		//���û�����ȼ����ײ ���Ǿ���Ҫ�ṩһ�������������һ�����Ƿ����������ڲ�

		//*/

		///*
		//���� ���ǻ�����ȷ�����������������ĵ�
		//ֻ����������Ҫ����洢һЩ��Ϣ
		//������Ǳ��������������Ʒ�˹������ε���������εĵ�

		//�����ĸ�����
		//�ֱ𱣴��������㵥���ε� A B������ԭ����ͼ���ϵ�λ��

		//���Ǿ�����������ȷ�����������������ĵ�
		//*/

		///* Convex Combination ͹���
		//S��һ��͹��
		//CH(S) = �Ʀ�iPi = ��1P1 + ... + ��nPn
		//where Pi��S , �ˡ�R
		//and �Ʀ�i = 1
		//where ��i >= 0
		//*/

		///*
		//����2D��˵
		//CH(S) = ��1P1 + ��2P2 Ҳ����˵�߶��ϵĵ�������߶ζ˵�����ض�ϵ������ʾ
		//���Q���㷨��ֹ�׶ε���������ԭ������ĵ�
		//��ô��Q��ԭ�������һ����ֱ��Q��λ�ڵ��߶�
		//L = A - B
		//Q * L = 0

		//��Q�ļ�Ȩ���Ӵ�����Եõ�
		//(��1*A + ��2 * B) * L = 0;
		//����
		//��1 + ��2 = 1
		//���������̿��Եõ� ��1 ��2

		//��ü�Ȩ����֮�� �����������Ʒ�˹������״�ϵĵ�����������
		//Acloset = ��1 * As1 + ��2 * Bs1  ���������������γ��߶�AB��ͼ��S1�϶�Ӧ�ĵ�
		//Bcloset = ��1 * As2 + ��2 * Bs2
		//*/

		///*
		//����һЩ����������Ҫ���

		//1.������ɷ�˹������״�ϵ�A���B����һ����
		//��ôL��������
		//��˵����ԭ�������㲻�������ɷ�˹����״�ı���
		//�������ϵ�һ������
		//����A��B��support������ͬһ��λ��
		//���ǿ��Է���A��B��
		//if(L.zero)
		//{
		//Acloset = A.s1;//A��һ��support�� ������������������� �ֱ���A�ϵ�һ���� �� B�ϵ�һ����
		//Bcloset = A.s2;
		//}

		//2.�ڶ���������Ǽ�Ȩ���ӻ���ָ�ֵ�����
		//���Ϊ�� ��ζ����һ�����ɷ�˹�����support��Ϊ�����
		//if(lambda1 < 0)
		//{
		//Acloest = B.s1;
		//Bcloest = B.s2;
		//}
		//else if(lambda2 < 0)
		//{
		//Acloset = A.s1;
		//Bcloset = A.s2;
		//}
		//*/

		////α��
		//d = ;//choose a direction
		//	 //�ڵ������а���������
		//simplex.add(support(A, B, d));
		//simplex.add(support(A, B, -d));
		////��ʼ����
		//while (true)
		//{
		//	//��õ���������������ԭ������ĵ�
		//	p = closepointtoorigin(simplex.a, simplex.b);
		//	//����Ƿ�Ϊ0
		//	if (p.zero())
		//	{
		//		//˵������Ӵ�
		//		//����Ӵ���ζ�� ����Ϊ0  ��
		//		return false;
		//	}

		//	// p to origin ���µ�����
		//	//�淶����Ϊ������Ҫ����ͶӰ����
		//	d = p.negate().normallize();
		//	//����������ϻ��һ���µ����Ʒ�˹�����
		//	c = support(A, B, d);
		//	//�»�õĵ���֮ǰ�������Ǹ������ ���Ǹ��ӽ�ԭ������?
		//	Real dc = c.dot(d);
		//	Real da = simpex.a.dot(d);
		//	//tolerance �� ���� <= 0
		//	if (dc - da < tolerance)
		//	{
		//		//������Ǹ��ӽ��� 
		//		//
		//		distance = dc;
		//		return true;
		//	}

		//	if (distance(a, ori) < distance(b, ori))
		//	{
		//		b = c;
		//	}
		//	else
		//	{
		//		a = c;
		//	}

		//}




		////----------------------------------��͸����봩͸����----------------------------------------

		///*
		//�������͹������Ʒ�˹������״����ԭ��
		//���������������ཻ��
		//���Ʒ�˹������״��ԭ�����̾��뼴Ϊ��͸���
		//ͬ����
		//������㵽ԭ���γɵ��������Ǵ�͸����
		//*/

		////EPA Expanding Polytope Algorigthm ���Ŷ������㷨
		///*
		//����Ҫ�����ɷ�˹������״���洴��һ�������� ������
		//������������ ֱ��ԭ�㵽�������������ı߾������ɷ�˹������״�ı�
		//���ŵķ�������
		//���ݶ����嵽ԭ�����̾��벻�����Ŷ����
		//ͨ������ ���ǲ����Ķ�����������ɷ�˹������״��Զ��������
		//�����Ϳ�����ô�͸��Ⱥʹ�͸����
		//EPA�ĵ����ο�������������
		//*/
		//void EPA()
		//{
		//	/*
		//	������
		//	EPA�㷨���ʼ��ʱ��
		//	��Ҫһ����ʼ���ĵ�����
		//	���ǿ�������GJK�㷨��ֹʱ��ĵ�������ΪEPA�㷨�ĳ�ʼ������

		//	EPA�㷨��Ҫ�����ĵ����� ��2d�������������
		//	��3d�������������
		//	*/
		//	simplex = ;//GJK�㷨��ֹ�ĵ�����
		//	while (true)
		//	{
		//		//��õ������Ͼ���ԭ������ı�
		//		Edge e = findcloseedge(s);
		//		//���һ���µ����ɷ�˹����� ��������ߵķ���
		//		vector p = support(A, B, e.normal);

		//		//������ڵ��������ǲ��Ǿ�����ԭ������ı�
		//		//�����ֵΪ0 ���� �ǳ�С�Ļ� ���Ƕ�������Ϊ�����������
		//		double d = p.dot(e.normal);
		//		if (d - e.distance < tolerance)
		//		{
		//			//tolerace �Ǹ��������� ex 1e-6
		//			//�������Ⱦ��Ȼ�ҪС
		//			//��ô���ǿ�����Ϊ���ǲ��ܰѵ�������չ�ĸ�Զ��
		//			//��ô���Ǿͻ���˽��
		//			normal = e.normal;
		//			depth = d;
		//		}
		//		else
		//		{
		//			//����û���ҵ�
		//			//�������Ǽ�����չ������
		//			//�������ӵ���Ӧ���±�λ�� 
		//			//simplex ���Ҫ�����������Ľṹ
		//			//ΪʲôҪ��ӵ����±���ص�λ����
		//			/*
		//			��Ϊ���ǲ�Ӧ���ƻ��ߵ�˳��
		//			����˵ ����������������֮��
		//			���ǻ����ܰ�˳������ȡ��
		//			*/
		//			simplex.insert(p, edge.index);
		//		}
		//	}
		//}

		////??�ⲻ���Ǳ���ô...
		////��tm���Ǳ�������С��
		///*

		//��������Ż�
		//�ò���ÿ�μ��㶼���¼������еı�
		//���Ѿ�������ıߵ�ֵ��������������
		//������ȱ�ǵ��ڴ�
		//*/
		//void find_close_edge()
		//{
		//	edge closet;
		//	//�ѱߵľ��������??
		//	closet.distance = Real_max_value;
		//	// s is the passed in simplex
		//	for (size_t i = 0; i < s.length; ++i)
		//	{
		//		//???
		//		size_t j = i + 1;
		//		j == s.length ? 0 : i + 1;//??
		//								  //get the current point and the next one
		//		vector a = s.get(i);
		//		vector b = s.get(j);
		//		//creat the edge vector
		//		vector e = b - a;
		//		//get vector OA
		//		vector oa = a - origin;
		//		//get the vector from edge towars the origin
		//		vector n = tripleProduct(e, oa, e);//����Զ��ԭ���Ǹ�����ķ�����
		//										   //normalize the vector
		//		n.mormalize();
		//		//calculate the distance from the origin to the dege
		//		// oooooooo  ����ԭ�㵽�����ߵľ���
		//		double d = n.dot(a);
		//		//check the distance against the other distance
		//		if (d < closet.distance)
		//		{
		//			closet.distance = d;//ԭ���������ߵľ���
		//			closet.normal = n;//����ߵķ�����
		//			closet.index = j;//�ڼ��������
		//		}
		//	}
		//	//return the closet edge we found
		//	return closet;
		//}




		///*
		//��Ҫע��ĵط�

		//�㷨��ʹ�����ػ��õ�һ�����ڱ���ԭ�㷽��ķ�����
		//���� ��С�Ļ��߽Ӵ���ײ������� ʹ�����ػ����ܻ�ʹEPA�㷨��������
		//����ԭ�������ߺܽ� ���ػ����ܷ���0����

		//�������ǿ����ñߵĴ�ֱ���������
		//if A = (x,y)
		//A.preproduct() = (-y,x) or (y,-x)
		//����������ϵ

		//���� ����ԭ��������߶�� ���Ƕ���ͨ����ֱ��ȡ�÷�����
		//if(clockwise) (y,-x);
		//else (-y,x);

		//���ֵ����ε���ת�������Ҫ
		//������ʹ�������Ч

		//���� ʲôʹ�����ε���ת����???
		//��ô����???

		//...

		//û�뵽
		//�������뵽����һ������
		//������Ҫִ��һ�������˷��Ϳ��Եó����
		//�����жϷ��� Ȼ������ ��ֱ��
		//*/


		///*
		//�Ľ� ???

		//��Ϊ�������ܴ� �ڴ�͸��Ƚ�Сʱ ����һ�㲻����EPA�㷨
		//����EPA�㷨һ��ʱGJK�㷨�Ĳ���
		//���ǻ���������ײ�����һ����(������״)�����м��;������
		//һ���õ���������֮��ľ��� ��ȥ������� �Ϳ��Եõ���������֮��ľ���
		//*

		////-------------------------------------------------------------------------------


} //namespace grid

#endif //_GRID_GEOMETRY_ALGORITHM_H_