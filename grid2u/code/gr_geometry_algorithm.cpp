//gr_geometry_algorithm.cpp

//TODO(zhong) : �ཻ��⺯���ǳ����� �������ǳ�������� ��Ӧ�ļ�������Ҳ�Ƚϲ�

//TODO(zhong) : ͬʱ����Ҳ��ҪһЩ��Լ򵥵Ĳ��Ժ���

#include"include/gr_geometry_algorithm.h"

namespace grid
{

	//------ֱ���߶������㷨----------

	void TwoPointsLineEquation(
		const Point& A, const Point& B,
		Real* coeA, Real* coeB, Real* coeC
	) {
		assert(coeA);
		assert(coeB);
		assert(coeC);
		*coeA = B.y - A.y;
		*coeB = A.x - B.x;
		*coeC = -(*coeA * A.x + *coeB * A.y);
	}

	void PointNormalLineEquation(
		const Point& P,
		const Vector2& normal,
		Real* coeA, Real* coeB, Real* coeC
	) {
		assert(coeA);
		assert(coeB);
		assert(coeC);
		*coeA = normal.x;
		*coeB = normal.y;
		*coeC = normal.x * P.x + normal.y * P.y;
	}

	bool SegmentStraddeleTest(
		const Segment& segment_1,
		const Segment& segment_2
	) {
		Vector2 AB(segment_1.end - segment_1.start);
		Vector2 AC(segment_2.start - segment_1.start);
		Vector2 AD(segment_2.end - segment_1.start);
		return AB.CrossProduct(AC) * AB.CrossProduct(AD) < 0.0;
	}

	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Point& LineC,
		const Point& LineD
	) {
		Vector2 CD = LineD - LineC;
		Vector2 CA = A - LineC;
		Vector2 CB = B - LineC;
		return CD.CrossProduct(CA) * CD.CrossProduct(CB) < 0.0;
	}

	//����������ֱ���ϵ����
	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Line& line
	) {
		Vector2 CD = VectorVerticalProduct(line.normal);
		Vector2 CA = A - line.origin;
		Vector2 CB = B - line.origin;
		return CD.CrossProduct(CA) * CD.CrossProduct(CB) < 0.0;
	}

	bool VectorRegionTest(
		const Vector2& OA, const Vector2& OB,
		const Vector2& OC
	) {
		Real rotation_direction = OA.CrossProduct(OB);
		return
			rotation_direction * OA.CrossProduct(OC) > 0.0 &&
			rotation_direction * OB.CrossProduct(OC) < 0.0;
	}

	//-----��Զ���㷨------

	int GetFartestPointIndex(
		const Vector2 direction,
		int size, const Point* point_array
	) {
		if (size <= 0)
			return -1;
		assert(point_array);
		int index = 0;
		Real distance = 0.0;
		Real max_distance = direction.DotProduct(point_array[index]);
		for (int i = 1; i < size; ++i)
		{
			distance = direction.DotProduct(point_array[i]);
			if (distance > max_distance) {
				index = i;
				max_distance = distance;
			}
		}
		return index;
	}

	Point GetFartestPoint(
		const Vector2 direction,
		int size, const Point* point_array
	) {
		int index = GetFartestPointIndex(direction, size, point_array);
		if (index != -1) {
			return point_array[index];
		}
		else {
			assert(index != -1);
			return Point(0, 0);
		}
	}

	PointNode* GetFartestPointPtr(
		const Vector2& direction,
		const PointSet& pointset
	) {
		if (pointset.Empty())
			return nullptr;
		PointNode* max_ptr = pointset.HeadPtr();
		Real distance = 0.0;
		Real max_distance = direction.DotProduct(max_ptr->point);
		PointNode* work = pointset.HeadPtr();
		do
		{
			distance = direction.DotProduct(work->point);
			if (distance > max_distance) {
				max_ptr = work;
				max_distance = distance;
			}
			work = work->next;
		} while (work != pointset.HeadPtr());
		return max_ptr;
	}

	Point GetFartestPoint(
		const Vector2& direction,
		const PointSet& pointset
	) {
		PointNode* work = GetFartestPointPtr(direction, pointset);
		if (work) {
			return work->point;
		}
		else {
			assert(work);
			return Point(0, 0);
		}
	}

	//-----�����㷨-------

	Real CalculatePointsDistance(
		const Point& point_1,
		const Point& point_2
	) {
		return (point_2 - point_1).Length();
	}

	Real CalculatePointToSegmentDistance(
		const Point& point,
		const Segment& segment
	) {
		if (segment.IsStandard()) {
			Vector2 PA = segment.start - point;
			Vector2 PB = segment.end - point;
			Vector2 normal = VectorVerticalProduct(segment.end - segment.start);
			Real distance = std::fabs(normal.DotProduct(Vector2(segment.start - point)));
			if (PA.DotProduct(normal) * PB.DotProduct(normal) <= 0.0) {
				return distance;
			}
			else {
				return std::sqrt(std::fmin(PA.LengthSquare(), PB.LengthSquare()));
			}
		}
		else {
			return (point - segment.start).Length();
		}
	}

	//----ͶӰ�����㷨-------

	Real CalculateProjectionLength(
		const Vector2& axis,
		int size, const Point* point_array,
		Real* min, Real* max,
		int* min_index, int* max_index
	) {
		assert(point_array);
		Real length = 0.0;
		int max_i = 0;
		int min_i = 0;
		Real max_length = 0.0;
		Real min_length = REAL_MAX;//REL_MAX
		for (int i = 0; i < size; ++i)
		{
			length = axis.DotProduct(point_array[i]);
			if (length > max_length) {
				max_length = length;
				max_i = i;
			}
			if (length < min_length) {
				min_length = length;
				min_i = i;
			}
		}
		if (min) *min = min_length;
		if (max) *max = max_length;
		if (min_index) *min_index = min_i;
		if (max_index) *max_index = max_i;
		return max_length - min_length;
	}

	Real CalculateProjectionLength(
		const Vector2& axis, const Polygon& polygon,
		Real* min, Real*max
	) {
		return CalculateProjectionLength(
			axis, polygon.n, polygon.vertex, min, max
		);
	}

	Real CalculateProjectionLength(
		const Vector2& axis, const Box& box,
		Real* min, Real* max
	) {
		Real min_length = 0.0;
		Real max_length = 0.0;
		min_length = max_length = box.center.DotProduct(axis);
		Real half = std::fabs(
			box.direction_x.DotProduct(axis)) +
			std::fabs(box.direction_y.DotProduct(axis));
		min_length -= half;
		max_length += half;
		if (min) *min = min_length;
		if (max)*max = max_length;
		return max_length - min_length;
	}

	Real CalculateProjectionLength(
		const Vector2& axis,
		const PointSet& pointset,
		Real* min, Real *max,
		PointNode* min_nptr, PointNode* max_nptr
	) {
		assert(!pointset.Empty());
		Real length = 0.0;
		PointNode* max_ptr = nullptr;
		PointNode* min_ptr = nullptr;
		Real max_length = 0.0;
		Real min_length = REAL_MAX;//REL_MAX
		PointNode* work = pointset.HeadPtr();
		do {
			length = axis.DotProduct(work->point);
			if (length > max_length) {
				max_length = length;
				max_ptr = work;
			}
			if (length < min_length) {
				min_length = length;
				min_ptr = work;
			}
			work = work->next;
		} while (work != pointset.HeadPtr());
		if (min) *min = min_length;
		if (max) *max = max_length;
		min_ptr = min_ptr;
		max_ptr = max_ptr;
		return max_length - min_length;
	}

	Real CalculateProjectionLength(
		const Vector2& axis,
		const Circle& circle,
		Real* min, Real* max
	) {
		if (min) *min =
			circle.center.DotProduct(axis) - circle.radius;
		if (max) *max =
			circle.center.DotProduct(axis) + circle.radius;
		return 2.0 * circle.radius;
	}


	Real CalculateProjectionLength(
		const Vector2& axis, const Triangle& triangle,
		Real* min, Real* max
	){
		return CalculateProjectionLength(
			axis, 3,triangle.vertex, min, max
		);
	}


	Real CalculateProjectionLength(
		const Vector2& axis, const Rectangle& rectangle,
		Real* min, Real* max
	) {
		Polygon polyrec(rectangle);
		return CalculateProjectionLength(
			axis, polyrec, min, max
		);
	}

	//--------���������㷨-------

	void FindHull(
		const PointSet& Sk,
		PointNode* P, PointNode* Q,
		PointSet* polygon)
	{
		assert(polygon);
		//�Ӽ���Sk�м���Ѱ�����͹���ĵ�
		//Sk�еĵ����߶�PQ�ұߵĵ�
		//���������û�е���˳�
		if (Sk.Empty()) {
			return;
		}
		else {
			//�Ӹ���������Ѱ�Ҿ����߶�PQ��Զ�ĵ�
			////�ҳ��õ�C
			//Point C(*Sk[max_distance_index]);
			PointNode* C(
				GetFartestPointPtr(
					VectorVerticalProduct(Q->point - P->point), Sk
				)
			);
			//���õ���ӵ�polygon�� λ����P���Q��֮��
			C = polygon->InsertAfter(P, C->point);
			//PQC���㼯�ֳ���������
			//S0λ��PQC֮��
			//S1��PC�Ҳ�
			//S2��CQ�Ҳ�
			//�ֱ��ҳ���Щ�� ��д���S1 S2
			PointSet S1;
			PointSet S2;
			Vector2 PC(C->point - P->point);
			Vector2 CQ(Q->point - C->point);
			PointNode* work = Sk.HeadPtr();
			do
			{
				if (PC.CrossProduct(work->point - P->point) > 0.0) {
					S1.Push(work->point);
				}
				if (CQ.CrossProduct(work->point - C->point) > 0.0) {
					S2.Push(work->point);
				}
				work = work->next;
			} while (work != Sk.HeadPtr());
			
			//������S1 S2��Ѱ�� �������͹���ĵ�
			FindHull(S1, P, C, polygon);
			FindHull(S2, C, Q, polygon);
		}
	}

	void QuickHull(
		int size, Point* point_array,
		PointSet* pointset
	) {
		assert(size > 2);
		assert(point_array);
		assert(pointset);
		pointset->Clear();

		//��ʱ����ȷ���㼯��������������		
		//�ҳ��㼯����Զ��������AB ������ҵ�������
		//��x��ͶӰ ֻ��ҪѰ��x�����������Сֵ
		Point A(GetFartestPoint(Vector2(-1, 0), size, point_array));
		Point B(GetFartestPoint(Vector2(1, 0), size, point_array));
		
		//��ӵ�͹����
		PointNode* Aptr = pointset->Push(A);
		PointNode* Bptr = pointset->Push(B);

		//�߶�AB��ʣ�µĵ��Ϊ��������	
		Vector2 AB = B - A;
		Real cross = 0.0;
		//���߶�AB�ұߵĵ�ΪS1
		//���߶�AB���ĵ�ΪS2
		PointSet S1;
		PointSet S2;
		for (int i = 0; i < size; ++i)
		{	//������	
			cross = AB.CrossProduct(point_array[i] - A);
			if (cross > 0.0) {//������>�� ��ӵ�S1
				S1.Push(point_array[i]);
			}
			else if (cross < 0.0) {//������<0 ��ӵ�S2
				S2.Push(point_array[i]);
			}
		}

		//��S1��S2�м���Ѱ�����͹���ĵ�
		FindHull(S1, Aptr, Bptr, pointset);
		FindHull(S2, Bptr, Aptr, pointset);

		//ִ�����֮�� polygon���ǿ���ʹ�õ�͹��
	}

	bool IsConvexPolygon(
		const Polygon& polygon
	) {
		if (polygon.n < 3)
			return false;
		Vector2 AB(polygon.vertex[0] - polygon.vertex[polygon.n - 1]);
		Vector2 BC(polygon.vertex[1] - polygon.vertex[0]);
		Real clock_wise = AB.CrossProduct(BC);
		if (clock_wise == 0.0) {
			return false;
		}
		for (int i = 1; i < polygon.n - 1; ++i) {
			AB = BC;
			BC = polygon.vertex[i + 1] - polygon.vertex[i];
			if (AB.CrossProduct(BC) * clock_wise <= 0.0) {
				return false;
			}
		}
		return true;
	}


	//��Χ���㷨
	Rectangle RectangleBounding(const Polygon& polygon)
	{
		Real x_min = polygon.vertex[0].x;
		Real x_max = polygon.vertex[0].y;
		Real y_min = polygon.vertex[0].x;
		Real y_max = polygon.vertex[0].y;
		for (int i = 1; i < polygon.n; ++i)
		{
			if (polygon.vertex[i].x < x_min)
				x_min = polygon.vertex[i].x;
			else if (polygon.vertex[i].x > x_max)
				x_max = polygon.vertex[i].x;
			if (polygon.vertex[i].y < y_min)
				y_min = polygon.vertex[i].y;
			else if (polygon.vertex[i].y > y_max)
				y_max = polygon.vertex[i].y;
		}
		return Rectangle(x_min, y_min, x_max - x_min, y_max - y_min);
	}

	//�ཻ����

	bool IsPointSegmentCollineation(
		const Point& point,
		const Segment& segment
	) {
		return (segment.start - point).CrossProduct(segment.end - point) == 0.0;
	}

	bool IsPointOnSegment(
		const Point& point,
		const Segment& segment
	) {
		if (IsPointSegmentCollineation(point, segment)) {
			if ((segment.start - point).DotProduct(segment.end - point) <= 0.0) {
				return true;
			}
			else {
				return false;
			}
		}
		else {
			return false;
		}
	}

	bool IsPointInTriangle(
		const Point& point,
		const Triangle& triangle
	) {
		Vector3 coefficient(1.0, point.x, point.y);;
		Matrix3 matrix(
			1.0, 1.0, 1.0,
			triangle.vertex[0].x, triangle.vertex[1].x, triangle.vertex[2].x,
			triangle.vertex[0].y, triangle.vertex[1].y, triangle.vertex[2].y
		);
		coefficient = matrix.Invert() * coefficient;
		if (
			coefficient.x > 0.0 && coefficient.x < 1.0 &&
			coefficient.y > 0.0 && coefficient.y < 1.0 &&
			coefficient.z > 0.0 && coefficient.z < 1.0
			) {
			return true;
		}
		else {
			return false;
		}
	}


	bool IsPointOnTriangle(
		const Point& point,
		const Triangle& triangle
	) {
		//�Ҿ���ֱ�ӱ��������߾Ϳ���
		for (int i = 0; i < 2; ++i){
			if (IsPointOnSegment(point, triangle.Side(i))) {
				return true;
			}
		}
		return false;
	}

	bool IsPointOutTriangle(
		const Point& point,
		const Triangle& triangle
	) {
		return !(
			IsPointInTriangle(point, triangle) || 
			IsPointOnTriangle(point, triangle)
			);
	}

	bool IsRectangleIntersection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	) {
		return
			rectangle_1.x < rectangle_2.x + rectangle_2.width &&
			rectangle_1.x + rectangle_1.width > rectangle_2.x &&
			rectangle_1.y < rectangle_2.y + rectangle_2.height&&
			rectangle_1.y + rectangle_1.height > rectangle_2.y;
	}



		//------------------------------//

	void PointToPointIntersectionDetection(
		const Point& point_1,
		const Point& point_2,
		IntersectionResult *result
	) {		
		assert(result);
		result->closet_point[0] = point_1;
		result->closet_point[1] = point_2;
		result->intersection_vector = (point_2 - point_1);
		result->intersection_depth = result->intersection_vector.Length();
		if (result->intersection_vector.IsZero()) {
			result->state = GeometryIntersectionState::INTERSECTION;
		}		
		else{
			result->intersection_vector.Normalize();
			result->state = GeometryIntersectionState::SEPERATION;		
		}
	}

	void PointToSegmentIntersectionDetection(
		const Point& point,
		const Segment& segment,
		IntersectionResult* result
	){
		assert(result);
		Vector2 PA = segment.start - point; //�������ϵ�һ��
		Vector2 PB = segment.end - point;	//�������ϵ�����һ��
		Vector2 AB = segment.end - segment.start;
		Vector2 normal = VectorTripleProduct(AB, PA, AB);
		if (normal.IsZero()) {
			//˵������
			if (PA.DotProduct(PB) > 0.0) {
				result->state = GeometryIntersectionState::SEPERATION;
				result->closet_point[0] = point;
				if (PA.Length() < PB.Length()) {
					//���˵��P����A���					
					result->closet_point[1] = segment.start;
					result->intersection_vector = PA.Normalize();
					result->intersection_depth = PA.Length();
				}
				else{//Ҳ����B��						
					result->closet_point[1] = segment.end;
					result->intersection_vector = PB.Normalize();
					result->intersection_depth = PB.Length();
				}
			}
			else {
				result->state = GeometryIntersectionState::INTERSECTION;
				result->closet_point[0] = result->closet_point[1] = point;
				result->intersection_depth = 0.0;
				result->intersection_vector.Set(0.0, 0.0);
			}
		}
		else{
			result->state = GeometryIntersectionState::SEPERATION;
			result->intersection_vector = normal.Normalize();
			result->intersection_depth = PA.DotProduct(
				result->intersection_vector
			);		
			result->closet_point[0] = point;
			result->closet_point[1] = 
				point + normal * result->intersection_depth;
		}		
	}

	void PointToTriangleIntersectionDetection(
		const Point& point,
		const Triangle& triangle,
		IntersectionResult* result
	) {
		//����һЩ����
		assert(result);
		//�������������β���һ����Ч�������� ��ô�㷨������
		//������Ҳ�����������Ϊ Ӧ�ñ���
		Vector3 coefficient(1.0, point.x, point.y);;
		Matrix3 matrix(
			1.0,1.0,1.0,
			triangle.vertex[0].x,triangle.vertex[1].x,triangle.vertex[2].x,
			triangle.vertex[0].y,triangle.vertex[1].y,triangle.vertex[2].y
		);
		coefficient = matrix.Invert() * coefficient;
		//�����Ҿ������ϵ��
		//������Щϵ�� �����ж����еĶ����
		if (
			coefficient.x > 0.0 && coefficient.x < 1.0 &&
			coefficient.y > 0.0 && coefficient.y < 1.0 &&
			coefficient.z > 0.0 && coefficient.z < 1.0
		){
			//�����������ڲ�
			result->state = GeometryIntersectionState::CONTAINT;
		}
		else if (
			IsZeroReal(coefficient.x) || 
			IsZeroReal(coefficient.y) || 
			IsZeroReal(coefficient.z)
			) {
			//������������
			result->state = GeometryIntersectionState::INTERSECTION;
			result->closet_point[0] = result->closet_point[1] = point;
			result->intersection_depth = 0.0;
			result->intersection_vector.Set(0.0, 0.0);
		}
		else {
			//������������
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}

	//����һ���ȽϷ����ʵ��
	//Ӧ�÷ǳ��� ȫ����ifelse ���Ǵ���Ƚ϶�
	void PointToRectangleIntersectionDetection(
		const Point& point,
		const Rectangle& rectangle,
		IntersectionResult* result
	) {
		assert(result);
		result->closet_point[0] = point;
		if (
			point.x > rectangle.x					&& 
			point.x < rectangle.x + rectangle.width &&
			point.y > rectangle.y					&& 
			point.y < rectangle.y + rectangle.height
		){
			result->state = GeometryIntersectionState::CONTAINT;
			result->closet_point[1].Set(
				point.x < rectangle.x + 0.5*rectangle.width ? 
				rectangle.x : rectangle.x + rectangle.width,
			point.y < rectangle.y + 0.5*rectangle.height ?
				rectangle.y : rectangle.y + rectangle.height
			);
			if (std::fabs(point.x - result->closet_point[1].x) <
				std::fabs(point.y - result->closet_point[1].y)) {
				result->closet_point[1].y = point.y;
			}
			else {
				result->closet_point[1].x = point.x;
			}
		}
		else {
			result->state = GeometryIntersectionState::SEPERATION;
			result->closet_point[1].Set(
				point.x <= rectangle.x ? rectangle.x : rectangle.x + rectangle.width,
				point.y <= rectangle.y ? rectangle.y : rectangle.y + rectangle.height
			);
		}
		result->intersection_depth =
			(result->closet_point[1] - result->closet_point[0]).Length();
		if (result->intersection_depth == 0.0) {
			result->state = GeometryIntersectionState::INTERSECTION;
			result->intersection_vector.Set(0.0, 0.0);
		}
		else {
			result->intersection_vector =
				(result->closet_point[1] - result->closet_point[0]).Normalize();
		}		
	}

	//������ת��
	//ʵ���ر�� ���ھ��ε���ײ�жϼ���
	void PointToBoxIntersectionDetection(
		const Point& point,
		const Box& box,
		IntersectionResult* result
	){
		assert(result);
		Matrix3 rotate;
		Rectangle box_rectangle;
		box.ToRectangle(&box_rectangle, &rotate);
		PointToRectangleIntersectionDetection(
			point,
			box_rectangle,
			result
		);
		//���е����Զ�Ҫ�������ĵ���ת��ȥ
		rotate.Invert();
		rotate.Transform(&(result->closet_point[0] -= box.center));
		result->closet_point[0] += box.center;
		rotate.Transform(&(result->closet_point[1] -= box.center));
		result->closet_point[1] += box.center;
		rotate.Transform(&(result->intersection_vector));
	}

	//void PointToCircleIntersectionDetection(
	//	const Point& point,
	//	const Circle& circle,
	//	IntersectionResult* result
	//) {
	//	assert(result);

	//	result->intersection_depth = _real(
	//		(point - circle.center).Length() - circle.radius
	//	);

	//	result->intersection_vector[0] = (circle.center - point);
	//	if (result->intersection_vector[0].IsZero()) {
	//		result->intersection_vector[0].Set(1.0, 0.0);
	//	}
	//	else {
	//		result->intersection_vector[0].Normalize();
	//	}
	//	result->intersection_vector[1] = -result->intersection_vector[0];

	//	result->closet_point[0] = point;
	//	result->closet_point[1] =
	//		point +
	//		result->intersection_vector[0] * result->intersection_depth;

	//	if (result->intersection_depth > 0.0) {
	//		//����Բ��		
	//		result->intersection = false;
	//		result->contact = false;
	//		result->contain[0] = false;
	//		result->contain[1] = false;
	//	}
	//	else if (result->intersection_depth < 0.0) {
	//		//����Բnei
	//		result->intersection = false;
	//		result->contact = false;
	//		result->contain[0] = false;
	//		result->contain[1] = true;
	//	}
	//	else {
	//		//����Բ��
	//		result->intersection = true;
	//		result->contact = true;
	//		result->contain[0] = false;
	//		result->contain[1] = false;
	//	}
	//}

	void PointToArcIntersectionDetection(
		const Point& point,
		const Arc& arc,
		IntersectionResult* result
	) {
		//�ǳ���
		//�����жϵ��Ƿ���Բ��
		if (_real((point - arc.center).Length() - arc.radius) == 0.0) {
			Vector2 OA = arc.start - arc.center;;
			Vector2 OB = arc.end - arc.center;
			Vector2 OP = point - arc.center;
			if ((OA.CrossProduct(OB)*OA.CrossProduct(OP) > 0.0 && OA.CrossProduct(OB) * OB.CrossProduct(OP) < 0.0) ||
				point == arc.start || point == arc.end) {
				//���˵�����ڻ���
				result->closet_point[0] = result->closet_point[1] = point;
				result->intersection_depth = 0.0;
				result->intersection_vector.Set(0.0, 0.0);
				result->state = GeometryIntersectionState::INTERSECTION;
			}
			else {
				//��Բ�� �����ڻ���
				//�� ��������� �������������˵�����һ��
				result->state = GeometryIntersectionState::SEPERATION;
			}
		}
		else {
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}

	void PointToSectorIntersectionDetection(
		const Point& point,
		const Sector& sector,
		IntersectionResult* result
	) {
		assert(result);
		Vector2 OA = sector.start - sector.center;
		Vector2 OB = sector.end - sector.center;
		Vector2 OP = point - sector.center;
		if (_real((point - sector.center).Length() - sector.radius) < 0.0 &&
			OA.CrossProduct(OB) * OA.CrossProduct(OP) > 0.0 &&
			OA.CrossProduct(OB) * OB.CrossProduct(OP) < 0.0) {
			result->state = GeometryIntersectionState::CONTAINT;
		}
		else {
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}

	void PointToPolygonIntersectionDetection(
		const Point& point,
		const Polygon& polygon,
		IntersectionResult* result
	) {
		assert(result);

	}

	void SegmentIntersectionDetection(
		const Segment& segment_1,
		const Segment& segment_2,
		IntersectionResult* result
	) {
		//�����ų�ʵ��
		Rectangle rectangle_1(segment_1);
		Rectangle rectangle_2(segment_2);
		if (IsRectangleIntersection(rectangle_1, rectangle_2)) {
			//�໥����ʵ��
			Vector2 OP = segment_1.end - segment_1.start;
			Vector2 OA = segment_2.start - segment_1.start;
			Vector2 OB = segment_2.end - segment_1.start;
			if (OA.CrossProduct(OP) * OB.CrossProduct(OP) < 0.0) {
				OP = segment_2.end - segment_2.start;
				OA = segment_1.start - segment_2.start;
				OB = segment_2.end - segment_2.start;
				if (OA.CrossProduct(OP) * OB.CrossProduct(OP) < 0.0) {
					result->state = GeometryIntersectionState::INTERSECTION;
				}
				else {
					result->state = GeometryIntersectionState::SEPERATION;
				}
			}
			else {
				result->state = GeometryIntersectionState::SEPERATION;
			}
		}
		else {
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}



	void SegmentRectangleIntersectionDetection(
		const Segment& segment,
		const Rectangle rectangle,
		IntersectionResult* result
	) {
		//���������ų�
		Rectangle segment_rectangle(segment);
		if (IsRectangleIntersection(segment_rectangle, rectangle)) {
			//����ʵ��
			//�жϾ����ĸ������Ƿ����߶�����
			Vector2 OP = segment.end - segment.start;
			Vector2 OA(Vector2(rectangle.x, rectangle.y) - segment.start);
			Vector2 OB(Vector2(rectangle.x + rectangle.width, rectangle.y) - segment.start);
			if (OP.CrossProduct(OA)*OP.CrossProduct(OB) < 0.0) {
				result->state = GeometryIntersectionState::INTERSECTION;
			}
			else {
				Vector2& OC = OB;
				OC = Vector2(rectangle.x + rectangle.width, rectangle.y + rectangle.height) - segment.start;
				if (OP.CrossProduct(OA)*OP.CrossProduct(OB) < 0.0) {
					result->state = GeometryIntersectionState::INTERSECTION;
				}
				else {
					Vector2& OD = OB;
					OD = Vector2(rectangle.x, rectangle.y + rectangle.height) - segment.start;
					if (OP.CrossProduct(OA)*OP.CrossProduct(OB) < 0.0) {
						result->state = GeometryIntersectionState::INTERSECTION;
					}
					else {
						result->state = GeometryIntersectionState::SEPERATION;
					}
				}
			}
		}
		else {
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}

	void RectangleCircleIntersectionDetection(
		const Rectangle& rectangle,
		const Circle& circle,
		IntersectionResult* result
	) {
		//���㵽���εļ��
		//��þ����Ͼ���Բ������ĵ�
		PointToRectangleIntersectionDetection(
			circle.center,
			rectangle,
			result
		);
		Vector2 closet_point(result->closet_point[1]);
		//����õ㵽Բ�ĵľ�����뾶�Ƚ�
		if (result->intersection_depth > circle.radius) {
			result->state = GeometryIntersectionState::SEPERATION;
		}
		else if (result->intersection_depth < circle.radius) {
			result->state = GeometryIntersectionState::INTERSECTION;
		}
		else {
			result->state = GeometryIntersectionState::INTER_CONTACT;
		}
	}

	void BoxCircleIntersectionDetection(
		const Box& box,
		const Circle& circle,
		IntersectionResult* result
	) {
		//�ǳ���
		//���Ȼ��box��rectangle
		Rectangle box_rectangle;
		Matrix3 rotate;
		Circle rotate_circle(circle);
		//�Ҳ�Ӧ�û�������
		box.ToRectangle(&box_rectangle, &rotate);
		rotate.Transform(&(rotate_circle.center -= box.center));
		rotate_circle.center += box.center;
		//������ͨ�ľ���Բ����ײ���
		RectangleCircleIntersectionDetection(
			box_rectangle, rotate_circle, result
		);
		//�ѽ����ת��ȥ
		rotate.Invert();
		rotate.Transform(&(result->closet_point[0] -= box.center));
		result->closet_point[0] += box.center;
		rotate.Transform(&(result->closet_point[1] -= box.center));
		result->closet_point[1] += box.center;
		rotate.Transform(&(result->intersection_vector));
	}

	//������������ͶӰ�ᷨ��
	//ûɶ�Ľ�
	//����Ҳ�ͽ���ʵ��һ�°�
	void RectangleBoxIntersectionDetection(
		const Rectangle& rectangle,
		const Box& box,
		IntersectionResult* result
	) {
		//���ǽ���Ҫȷ���ĸ���

	}

	//OBB
	void BoxIntersectionDetection(
		const Box& box_1,
		const Box& box_2,
		IntersectionResult* result
	) {
		Real box_1_proj_min = 0.0;
		Real box_1_proj_max = 0.0;
		Real box_2_proj_min = 0.0;
		Real box_2_proj_max = 0.0;
		Vector2 *axis[4] = {
			&((const_cast<Box&>(box_1)).direction_x),
			&((const_cast<Box&>(box_1)).direction_y),
			&((const_cast<Box&>(box_2)).direction_x),
			&((const_cast<Box&>(box_2)).direction_y)
		};
		int count = 4;
		if (box_1.direction_x.CrossProduct(box_2.direction_x) == 0.0 ||
			box_1.direction_x.DotProduct(box_2.direction_x) == 0.0) {
			count = 2;
		}			
		for (int i = 0; i < count; ++i) {
			CalculateProjectionLength(*axis[i], box_1, &box_1_proj_min, &box_1_proj_max);
			CalculateProjectionLength(*axis[i], box_2, &box_2_proj_min, &box_2_proj_max);
			if (box_1_proj_min <= box_2_proj_max && box_2_proj_min <= box_1_proj_max) {
				continue;
			}
			else {
				result->state = GeometryIntersectionState::SEPERATION;
				return;
			}
		}
		result->state = GeometryIntersectionState::INTERSECTION;		
	}


	void CircleIntersectionDetection(
		const Circle& circle_1,
		const Circle& circle_2,
		IntersectionResult* result
	) {
		Real center_distance = (circle_1.center - circle_2.center).Length();
		if (center_distance <= 0.0 &&center_distance < std::fabs(circle_1.radius - circle_2.radius)) {//�ں�
			result->state = GeometryIntersectionState::CONTAINT;
		}
		else if (_real(center_distance - std::fabs(circle_1.radius - circle_2.radius)) == 0.0) {//����
			result->state = GeometryIntersectionState::INTER_CONTACT;
		}
		else if (center_distance > std::fabs(circle_1.radius - circle_2.radius) && center_distance < circle_1.radius + circle_2.radius) { //�ཻ
			result->state = GeometryIntersectionState::SEPERATION;
		}
		else if (_real(center_distance - (circle_1.radius + circle_2.radius)) == 0.0) { //����
			result->state = GeometryIntersectionState::OUTER_CONTACT;
		}
		else { //����
			result->state = GeometryIntersectionState::SEPERATION;
		}
	}

	void SegmentCircleIntersectionDetection(
		const Segment& segment,
		const Circle& circle,
		IntersectionResult* result
	) {
		//���ȼ���Բ�ĵ��߶εľ���
		Real distance = CalculatePointToSegmentDistance(circle.center, segment);
		//��뾶�Ƚ�
		if (distance > circle.radius) {
			result->state = GeometryIntersectionState::SEPERATION;
		}
		else if (distance == circle.radius) {
			//�Ӵ�
			result->state = GeometryIntersectionState::OUTER_CONTACT;
		}
		else {
			//�߶ε�һͷ��Բ���ڲ� 
			//���Ǵ�ʱ��Ҫȷ������һͷ��Բ���ڲ������ⲿ
			//����ԭ�㵽��Զ��ľ���
			//distance = CalculatePointToSegmentMaxDistance(
			//	circle.center, segment
			//);
			if (distance < circle.radius) {
				//�߶���Բ��
				result->state = GeometryIntersectionState::CONTAINT;
			}
			else if (distance == circle.radius) {
				//�߶���Բ�ڽӴ�
				result->state = GeometryIntersectionState::INTER_CONTACT;
			}
			else {
				//�߶���ԭ�ཻ
				result->state = GeometryIntersectionState::INTERSECTION;
			}
		}
	}





	void CircleSectorIntersectionDetection(
		const Circle& circle,
		const Sector& sector,
		IntersectionResult* result
	) {
		//�����������ų�
		//�ж�Բ���Ƿ�������������
		PointToSectorIntersectionDetection(
			circle.center,
			sector,
			result
		);
		//����λ�ù�ϵ
		//�����Բ��
		if (result->state == GeometryIntersectionState::CONTAINT) {
			//���ǱȽϸ���
		}
		else {

		}



	}



	Vector2 GetClosetPointOnPolygon(
		const Point& point,
		const Polygon& polygon
	) {
		Vector2 closet(polygon.vertex[0]);
		Real closet_distance = (point - polygon.vertex[0]).LengthSquare();
		Real distance = 0.0;
		for (int i = 1; i < polygon.n; ++i) {
			distance = (point - polygon.vertex[i]).LengthSquare();
			if (distance < closet_distance) {
				closet_distance = distance;
				closet = polygon.vertex[i];
			}
		}
		return closet;
	}



	bool CirclePolygonSAT(
		const Circle& circle,
		const Polygon& polygon
	) {
		//������һ����
		//����Բ��һ������ �����������Ȼ����Բ�йص�������
		Real poly_proj_min = 0.0;
		Real poly_proj_max = 0.0;
		Real cir_proj_min = 0.0;
		Real cir_proj_max = 0.0;
		Vector2 axis(0.0, 0.0);
		//������Ҫ��ö�������ж�������Բ������ĵ�
		Vector2 closet = GetClosetPointOnPolygon(
			circle.center,
			polygon
		);
		axis = circle.center - closet;
		if (axis.IsZero()) {
			return false;
		}
		CalculateProjectionLength(axis , circle, &cir_proj_min, &cir_proj_max);
		CalculateProjectionLength(axis, polygon, &poly_proj_min, &poly_proj_max);
		//����Ƿ��ص�
		if (poly_proj_min <= cir_proj_max && cir_proj_min <= poly_proj_max) {
			;
		}
		else {
			return true;
		}
		//Բ�μ�����֮��ʼ�������
		for (int i = 0; i < polygon.n; ++i) {
			axis = polygon.SideNormal(i);
			CalculateProjectionLength(axis,circle, &cir_proj_min, &cir_proj_max);
			CalculateProjectionLength(axis, polygon, &poly_proj_min, &poly_proj_max);
			//����Ƿ��ص�
			if (poly_proj_min <= cir_proj_max && cir_proj_min <= poly_proj_max) {
				continue;
			}
			else {
				return true;
			}
		}
		return false;
	}


} //namespace grid