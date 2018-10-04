//@header : gr_intersection.cpp
//@author : zhong

#include"include/gr_intersection.h"

namespace grid
{
	//PointIntersection

	bool IntersectionTest(
		const Point& point_1,
		const Point& point_2,
		bool REGION
	) {
		return point_1 == point_2;
	}

	bool IntersectionTest(
		const Point& point,
		const Line& line,
		bool REGION
	) {
		return IsVectorVertical(
			line.normal, point - line.origin
		);
	}

	bool IntersectionTest(
		const Point& point,
		const Ray& ray,
		bool REGION
	) {
		Line line(ray.origin, VectorVerticalProduct(ray.direction));
		if (IntersectionTest(point, line)) {
			Real projection = (point - ray.origin).DotProduct(ray.direction);
			return projection >= 0.0;
		}
		return true;
	}

	bool IntersectionTest(
		const Point& point,
		const Segment& segment,
		bool REGION
	) {
		Vector2 PA(segment.start - point);
		Vector2 PB(segment.end - point);
		return IsVectorParallel(PA, PB) && (PA).DotProduct(PB) <= 0.0;
	}

	bool IntersectionTest(
		const Point& point,
		const Triangle& triangle,
		bool REGION
	) {
		//Vector3 coefficient(1.0, point.x, point.y);;
		//Matrix3 matrix(
		//	1.0, 1.0, 1.0,
		//	triangle.vertex[0].x, triangle.vertex[1].x, triangle.vertex[2].x,
		//	triangle.vertex[0].y, triangle.vertex[1].y, triangle.vertex[2].y
		//);
		//coefficient = matrix.Invert() * coefficient;
		//return
		//	coefficient.x > 0.0 && coefficient.x < 1.0 &&
		//	coefficient.y > 0.0 && coefficient.y < 1.0 &&
		//	coefficient.z > 0.0 && coefficient.z < 1.0;
		//一种更为简单的测试方法
		//依次计算叉积
		Vector2 AB(triangle.vertex[1] - triangle.vertex[0]);
		Vector2 AP(point - triangle.vertex[0]);
		Real r = AB.CrossProduct(AP);
		AB = triangle.vertex[2] - triangle.vertex[1];
		AP = point - triangle.vertex[1];
		if (r * AB.CrossProduct(AP) < 0.0) {
			return false;
		}
		AB = triangle.vertex[0] - triangle.vertex[2];
		AP = point - triangle.vertex[2];
		if (r * AB.CrossProduct(AP) < 0.0) {
			return false;
		}
		return true;
	}

	bool IntersectionTest(
		const Point& point,
		const Rectangle& rectangle,
		bool REGION
	) {
		if (REGION) {
			return
				point.x >= rectangle.x &&
				point.x <= rectangle.x + rectangle.width &&
				point.y >= rectangle.y &&
				point.y <= rectangle.y + rectangle.height;
		}
		else {
			return
				(point.x == rectangle.x && point.y >= rectangle.y && point.y <= rectangle.y + rectangle.height) ||
				(point.x == rectangle.x + rectangle.width && point.y >= rectangle.y && point.y <= rectangle.y + rectangle.height) ||
				(point.y == rectangle.y && point.x >= rectangle.x && point.x <= rectangle.x + rectangle.width) ||
				(point.y == rectangle.y + rectangle.height && point.x >= rectangle.x && point.x <= rectangle.x + rectangle.width);
		}
	}

	bool IntersectionTest(
		const Point& point,
		const Box& box,
		bool REGION
	) {
		Rectangle rectangle;
		Matrix3 rotate;
		Point rotated_point(point);
		box.ToRectangle(&rectangle, &rotate);
		rotate.Transform(&rotated_point);
		return IntersectionTest(rotated_point, rectangle, REGION);
	}

	bool IntersectionTest(
		const Point& point,
		const Polygon& polygon,
		bool REGION
	) {
		assert(polygon.n > 2);
		int low = 0;
		int high = polygon.n;
		int mid = 0;
		bool cw = (polygon.vertex[1] - polygon.vertex[0]).CrossProduct(polygon.vertex[2] - polygon.vertex[1]) > 0.0;
		Triangle triangle;
		do {
			mid = (low + high) / 2;
			if (triangle.Set(polygon.vertex[0],polygon.vertex[mid],point).ClockWise() && cw) {
				low = mid;
			}
			else {		
				low = high;
			}
		} while (low + 1 < high);
		if (low == 0 || high == 0) {
			return false;
		}
		return triangle.Set(polygon.vertex[low], polygon.vertex[high], point).ClockWise() && cw;
	}

	bool IntersectionTest(
		const Point& point,
		const Arc& arc,
		bool REGION
	) {
		return
			IntersectionTest(point, 
				Circle(arc.center, arc.radius), 
				false) &&
			VectorRegionTest(arc.start - arc.center, 
				arc.end - arc.center, 
				point - arc.center);
	}

	bool IntersectionTest(
		const Point& point,
		const Sector& sector,
		bool REGION
	) {
		return
			IntersectionTest(point, 
				Circle(sector.center, sector.radius), 
				true) &&
			VectorRegionTest(sector.start - sector.center, 
				sector.end - sector.center, 
				point - sector.center);
	}

	bool IntersectionTest(
		const Point& point,
		const Circle& circle,
		bool REGION
	) {
		Real distance_square = (point - circle.center).LengthSquare();
		if (REGION) {
			return distance_square <= circle.radius * circle.radius;
		}
		else {
			return IsZeroReal(distance_square - circle.radius * circle.radius);
		}
	}

	//LineIntersection

	bool IntersectionTest(
		const Line& line,
		const Point& point,
		bool REGION
	) {
		return IntersectionTest(point, line, REGION);
	}

	bool IntersectionTest(
		const Line& line_1,
		const Line& line_2,
		bool REGION
	) {
		return !IsVectorParallel(line_1.normal, line_2.normal);
	}

	bool IntersectionTest(
		const Line& line,
		const Ray& ray,
		bool REGION
	) {
		Line line_ray(ray);
		Point point;
		if (IntersectionTest(line, line_ray)){
			line.IntersectionPoint(line_ray, &point);
			return (point - ray.origin).DotProduct(ray.direction) >= 0.0;
		}
		return false;
	}

	bool IntersectionTest(
		const Line& line,
		const Segment& segment,
		bool REGION
	) {
		return LineStraddleTest(segment.start, segment.end, line);
	}

	bool IntersectionTest(
		const Line& line,
		const Triangle& triangle,
		bool REGION
	) {
		Vector2 direction(VectorVerticalProduct(line.normal));		
		Real a = (triangle.vertex[0] - line.origin).DotProduct(direction);
		Real b = (triangle.vertex[1] - line.origin).DotProduct(direction);
		Real c = (triangle.vertex[2] - line.origin).DotProduct(direction);
		return a * b <= 0.0 || a * c <= 0.0 || b * c <= 0.0;
	}

	bool IntersectionTest(
		const Line& line,
		const Rectangle& rectangle,
		bool REGION
	) {
		return true;
	}

	bool IntersectionTest(
		const Line& line,
		const Box& box,
		bool REGION
	) {
		return true;
	}

	bool IntersectionTest(
		const Line& line,
		const Polygon& polygon,
		bool REGION
	) {
		Vector2 direction(VectorVerticalProduct(line.normal));
		Real r = (polygon.vertex[0] - line.origin).DotProduct(direction);
		for (int i = 1; i < polygon.n; ++i)
		{
			if ((polygon.vertex[i] - line.origin).DotProduct(direction) * r <= 0.0){
				return true;
			}
		}
		return true;
	}

	bool IntersectionTest(
		const Line& line,
		const Arc& arc,
		bool REGION
	) {
		return true;
	}

	bool IntersectionTest(
		const Line& line,
		const Sector& sector,
		bool REGION
	) {
		return true;
	}

	bool IntersectionTest(
		const Line& line,
		const Circle& circle,
		bool REGION
	) {
		return GrFabs((circle.center - line.origin).DotProduct(line.normal)) <= circle.radius;
	}

	//RayIntersectionTest

	bool IntersectionTest(
		const Ray& ray_1,
		const Ray& ray_2
	) {
		if (IsVectorParallel(ray_1.direction, ray_1.direction)) {
			return false;
		}
		return VectorRegionTest(
			ray_2.origin - ray_2.origin,
			ray_1.direction,
			ray_2.direction
		);
	}

	//SegmentIntersection

	bool IntersectionTest(
		const Segment& segment_1,
		const Segment& segment_2
	) {
		Rectangle rectangle_1(segment_1);
		Rectangle rectangle_2(segment_2);
		if (IntersectionTest(rectangle_1, rectangle_2)) {
			return SegmentStraddeleTest(segment_1, segment_2) &&
				SegmentStraddeleTest(segment_2, segment_1);
		}
		else {
			return false;
		}
	}

	//TriangleIntersection

	bool IntersectionTest(
		const Triangle& triangle_1,
		const Triangle& triangle_2
	) {
		Polygon polytri_1(triangle_1);
		Polygon polytri_2(triangle_2);
		return IntersectionTest(polytri_1, polytri_2);
	}

	//RectangleIntersection

	bool IntersectionTest(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	) {
		return
			rectangle_1.x < rectangle_2.x + rectangle_2.width  &&
			rectangle_1.x + rectangle_1.width > rectangle_2.x  &&
			rectangle_1.y < rectangle_2.y + rectangle_2.height &&
			rectangle_1.y + rectangle_1.height > rectangle_2.y ;
	}

	//BoxIntersection

	bool IntersectionTest(
		const Box& box_1,
		const Box& box_2
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
				return false;
			}
		}
		return true;
	}

	//PolygonIntersection

	bool IntersectionTest(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	) {
		Real poly_1_proj_min = 0.0;
		Real poly_1_proj_max = 0.0;
		Real poly_2_proj_min = 0.0;
		Real poly_2_proj_max = 0.0;
		Vector2 axis(0.0, 0.0);
		for (int i = 0; i < polygon_1.n + polygon_2.n; ++i) {
			if (i < polygon_1.n) {
				axis = polygon_1.SideNormal(i);
			}
			else {
				axis = polygon_2.SideNormal(i - polygon_1.n);
			}
			CalculateProjectionLength(axis, polygon_1.n, polygon_1.vertex, &poly_1_proj_min, &poly_1_proj_max);
			CalculateProjectionLength(axis, polygon_2.n, polygon_2.vertex, &poly_2_proj_min, &poly_2_proj_max);
			if (poly_1_proj_min <= poly_2_proj_max && poly_2_proj_min <= poly_1_proj_max) {
				continue;
			}
			else {
				return true;
			}
		}
		return false;
	}

	//ArcIntersection

	bool IntersectionTest(
		const Arc& arc_1,
		const Arc& arc_2
	) {
		Circle circle_1(arc_1.center, arc_1.radius);
		Circle circle_2(arc_2.center, arc_2.radius);
		if (IntersectionTest(circle_1, circle_2)) {
			Segment segment_1(arc_1.start, arc_2.end);
			Segment segment_2(arc_1.end, arc_2.start);
			if (IntersectionTest(segment_1, segment_2)) {
				return true;
			}
			else {
				segment_1.Set(arc_1.start, arc_2.start);
				segment_2.Set(arc_1.end, arc_2.end);
				if (IntersectionTest(segment_1, segment_2)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
		else {
			return false;
		}
	}

	//SectorIntersection

	bool IntersectionTest(
		const Sector& sector_1,
		const Sector& sector_2
	) {
		//首先实现圆与扇形的相交测试
		//然后依次判断是否互相相交
		//如果是 则相交
		//如果不是 则不相交
		return true;
	}

	//CircleIntersection

	bool IntersectionTest(
		const Circle& circle_1,
		const Circle& circle_2
	) {
		Real center_distance_square = 
			(circle_1.center - circle_2.center).LengthSquare();
		Real radius_difference_square =
			(circle_1.radius - circle_2.radius) * (circle_1.radius - circle_2.radius);
		Real radius_sum_square = 
			(circle_1.radius + circle_2.radius) * (circle_1.radius + circle_2.radius);
		return	center_distance_square >= radius_difference_square &&
			center_distance_square <= radius_sum_square;
	}

};