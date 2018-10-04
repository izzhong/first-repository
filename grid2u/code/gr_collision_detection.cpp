//@source : gr_collision_detection.cpp

#include"include/gr_geometry_algorithm.h"

namespace grid
{
	bool CircleSeparationDetection(
		const Circle& circle_1,
		const Circle& circle_2
	) {
		return
			(circle_1.center - circle_2.center).LengthSquare() >
			(circle_1.radius + circle_2.radius) *
			(circle_1.radius + circle_2.radius);
	}

	bool RectangleSeparationDetection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	) {
		return !(
			rectangle_1.x < rectangle_2.x + rectangle_2.width &&
			rectangle_1.x + rectangle_1.width > rectangle_2.x &&
			rectangle_1.y < rectangle_2.y + rectangle_2.height&&
			rectangle_1.y + rectangle_1.height > rectangle_2.y
			);
	}

	bool BoxSeparationDetection(
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
				return true;
			}
		}
		return false;
	}

	bool PolygonSeparationDetection(
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
				return false;
			}
		}
		return true;
	}


} //namespace grid