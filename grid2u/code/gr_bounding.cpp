//@source	:	gr_bounding.cpp
//@date		:	2018/9
//@author	:	zhong

#include"include/gr_bounding.h"

namespace grid
{
	void CenterBoundingCircle::CreateBounding(
		int size, const Point* point_array
	) {
		CreateCenterBoundingCircle(size, point_array, &circle_);
	}

	void CenterBoundingCircle::CreateBounding(
		const PointSet& pointset
	) {
		PointArray pointarry(pointset);
		CreateCenterBoundingCircle(
			pointarry.Size(), pointarry.Array(),
			&circle_
		);
	}

	void CenterBoundingCircle::CreateBounding(
		const Circle& circle
	) {
		circle_ = circle;
	}

	void CenterBoundingCircle::CreateBounding(
		const Triangle& triangle
	) {
		Point A(triangle.vertex[0]);
		Point B(triangle.vertex[1]);
		Point C(triangle.vertex[2]);
		Matrix2 coe;
		Vector2 c;
		PointNormalLineEquation(
			(A + B)*0.5, B - A,
			&coe.m[0], &coe.m[1], &c.x
		);
		PointNormalLineEquation(
			(A + C)*0.5, C - A,
			&coe.m[2], &coe.m[3], &c.y
		);
		if (coe.IsInvertible()) {
			circle_.center = coe.Invert()*c;
			circle_.radius = (A - circle_.center).Length();
		}
	}

	void CenterBoundingCircle::CreateBounding(
		const Rectangle& rectangle
	) {
		circle_.center = rectangle.Center();
		circle_.radius = 0.5 *  GrSqrt(
			rectangle.width * rectangle.width +
			rectangle.height * rectangle.height
		);
	}

	void CenterBoundingCircle::CreateBounding(
		const Box& box
	) {
		circle_.center = box.center;
		circle_.radius = GrSqrt(
			box.half_width * box.half_width +
			box.half_height * box.half_height
		);
	}

	void CenterBoundingCircle::CreateBounding(
		const Polygon& polygon
	) {
		CreateBounding(polygon.n, polygon.vertex);
	}

	void CenterBoundingCircle::CreateBounding(
		const Arc& arc
	) {
		circle_.center = (arc.start + arc.end) * 0.5;
		circle_.radius = (arc.start - circle_.center).Length();
	}

	void CenterBoundingCircle::CreateBounding(
		const Sector& sector
	) {
		circle_.center = sector.center;
		circle_.radius = sector.radius;
	}

	void CenterBoundingCircle::CreateCenterBoundingCircle(
		int size, const Point* point_array,
		Circle* cb
	) {
		assert(size > 0);
		assert(point_array);
		assert(cb);
		AABB aabb;
		aabb.CreateBounding(size, point_array);
		Rectangle rec(aabb.GetGeometry());
		Point center(rec.Center());
		Real radius = (Point(rec.x, rec.y) - center).Length();
		Circle circle(center, radius);
		Real radius_square = radius * radius;
		Real length_square = 0.0;
		for (int i = 0; i < size; ++i)
		{
			length_square = (point_array[i] - center).LengthSquare();
			if (length_square > radius_square)
				radius_square = length_square;
		}
		circle.radius = GrSqrt(radius_square);
		cb->center = center;
		cb->radius = circle.radius;
	}

	CenterBoundingCircle CombineBounding(
		const CenterBoundingCircle& cbc_1,
		const CenterBoundingCircle& cbc_2
	) {	
		Circle cir1 = cbc_1.GetGeometry();
		Circle cir2 = cbc_2.GetGeometry();
		Real center_distance = (cir1.center - cir2.center).Length();
		if (center_distance <= GrFabs(cir1.radius - cir2.radius)) {
			if (cir1.radius > cir2.radius) {
				return CenterBoundingCircle(cir1);
			}
			else {
				return CenterBoundingCircle(cir2);
			}
		}
		else{
			Circle circle;
			circle.radius = (center_distance + cir1.radius + cir2.radius) * 0.5;
			circle.center = cir1.center + (cir2.center - cir1.center) * 
				(circle.radius - cir1.radius) * (1.0 / center_distance);		
			return CenterBoundingCircle(circle);
		}
	}

	void AxisAlignedBoundingBox::CreateBounding(
		int size, const Point* point_array
	) {
		assert(size > 0);
		assert(point_array);
		Real x_min = point_array[0].x;
		Real x_max = point_array[0].y;
		Real y_min = point_array[0].x;
		Real y_max = point_array[0].y;
		for (int i = 1; i < size; ++i)
		{
			if (point_array[i].x < x_min)
				x_min = point_array[i].x;
			else if (point_array[i].x > x_max)
				x_max = point_array[i].x;
			if (point_array[i].y < y_min)
				y_min = point_array[i].y;
			else if (point_array[i].y > y_max)
				y_max = point_array[i].y;
		}
		rect_.x = x_min;
		rect_.y = y_min;
		rect_.width = x_max - x_min;
		rect_.height = y_max - y_min;
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const PointSet& point_set
	) {
		PointNode* work = point_set.HeadPtr();
		Real x_min = work->point.x;
		Real x_max = work->point.y;
		Real y_min = work->point.x;
		Real y_max = work->point.y;
		do
		{
			if (work->point.x < x_min)
				x_min = work->point.x;
			else if (work->point.x > x_max)
				x_max = work->point.x;
			if (work->point.y < y_min)
				y_min = work->point.y;
			else if (work->point.y > y_max)
				y_max = work->point.y;
			work = work->next;
		} while (work != point_set.HeadPtr());
		rect_.x = x_min;
		rect_.y = y_min;
		rect_.width = x_max - x_min;
		rect_.height = y_max - y_min;
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Circle& circle
	) {
		rect_.x = circle.center.x - circle.radius;
		rect_.y = circle.center.y - circle.radius;
		rect_.width = 2.0 * circle.radius;
		rect_.height = 2.0 * circle.radius;
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Triangle& triangle
	) {
		CreateBounding(3, triangle.vertex);
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Rectangle& rectangle
	) {
		rect_ = rectangle;
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Box& box
	) {
		Point vertex[4] =
		{
			box.center - box.direction_x * box.half_width -
			box.direction_y * box.half_height,
			box.center + box.direction_x * box.half_width -
			box.direction_y * box.half_height,
			box.center + box.direction_x * box.half_width +
			box.direction_y * box.half_height,
			box.center - box.direction_x * box.half_width +
			box.direction_y * box.half_height,
		};
		CreateBounding(4, vertex);
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Polygon& polygon
	) {
		CreateBounding(polygon.n, polygon.vertex);
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Arc& arc
	) {
		Point O(arc.center);
		Point A(O + Point(-arc.radius, 0.0));
		Point B(O + Point(0.0, -arc.radius));
		Point C(O + Point(arc.radius, 0.0));
		Point D(O + Point(0.0, arc.radius));
		Real x_max = 0.0;
		Real x_min = 0.0;
		Real y_min = 0.0;
		Real y_max = 0.0;
		if (LineStraddleTest(arc.start, arc.end, A, C)) {
			if (LineStraddleTest(A, arc.end, O, arc.start)) {
				x_max = C.x;
				x_min = GrMin(arc.start.x, arc.end.x);
			}
			else {
				x_min = A.x;
				x_max = GrMax(arc.start.x, arc.end.x);
			}
		}
		else {
			x_min = GrMin(arc.start.x, arc.end.x);
			x_max = GrMax(arc.start.x, arc.end.x);
		}
		if (LineStraddleTest(arc.start, arc.end, B, D)) {
			if (LineStraddleTest(B, arc.end, O, arc.start)) {
				y_max = D.x;
				y_min = GrMin(arc.start.y, arc.end.y);
			}
			else {
				y_min = B.y;
				y_max = GrMax(arc.start.y, arc.end.y);
			}
		}
		else {
			y_min = GrMin(arc.start.y, arc.end.y);
			y_max = GrMax(arc.start.y, arc.end.y);
		}
		rect_.x = x_min;
		rect_.y = y_min;
		rect_.width = x_max - x_min;
		rect_.height = y_max - y_min;
	}

	void AxisAlignedBoundingBox::CreateBounding(
		const Sector& sector
	) {
		Point O(sector.center);
		Point A(O + Point(-sector.radius, 0.0));
		Point B(O + Point(0.0, -sector.radius));
		Point C(O + Point(sector.radius, 0.0));
		Point D(O + Point(0.0, sector.radius));
		Real x_max = 0.0;
		Real x_min = 0.0;
		Real y_min = 0.0;
		Real y_max = 0.0;
		if (LineStraddleTest(sector.start, sector.end, A, C)) {
			if (LineStraddleTest(A, sector.end, O, sector.start)) {
				x_max = C.x;
				x_min = GrMin(GrMin(sector.start.x, sector.end.x),
					sector.center.x);
			}
			else {
				x_min = A.x;
				x_max = GrMax(GrMax(sector.start.x, sector.end.x),
					sector.center.x);
			}
		}
		else {
			x_min = GrMin(GrMin(sector.start.x, sector.end.x),
				sector.center.x);
			x_max = GrMax(GrMax(sector.start.x, sector.end.x),
				sector.center.x);
		}
		if (LineStraddleTest(sector.start, sector.end, B, D)) {
			if (LineStraddleTest(B, sector.end, O, sector.start)) {
				y_max = D.x;
				y_min = GrMin(GrMin(sector.start.x, sector.end.x),
					sector.center.x);
			}
			else {
				y_min = B.y;
				y_max = GrMax(GrMax(sector.start.x, sector.end.x),
					sector.center.x);
			}
		}
		else {
			y_min = GrMin(GrMin(sector.start.y, sector.end.y),
				sector.center.y);
			y_max = GrMax(GrMax(sector.start.y, sector.end.y),
				sector.center.y);
		}
		rect_.x = x_min;
		rect_.y = y_min;
		rect_.width = x_max - x_min;
		rect_.height = y_max - y_min;
	}

	AxisAlignedBoundingBox CombineBounding(
		const AxisAlignedBoundingBox& aabb_1,
		const AxisAlignedBoundingBox& aabb_2
	) {
		Real xmin = GrMin(aabb_1.rect_.x, aabb_2.rect_.x);
		Real xmax = GrMax(
			aabb_1.rect_.x + aabb_1.rect_.width, 
			aabb_2.rect_.x + aabb_2.rect_.width
		);
		Real ymin = GrMin(aabb_1.rect_.y, aabb_2.rect_.y);
		Real ymax = GrMax(
			aabb_1.rect_.y + aabb_1.rect_.height,
			aabb_2.rect_.y + aabb_2.rect_.height
		);
		Rectangle rectangle(xmin, ymin, xmax - xmin, ymax - ymin);
		return AxisAlignedBoundingBox(rectangle);
	}

	void OrientedBoundingBox::CreateBounding(
		int size, const Point* pointarray
	) {
		Vector2 directionx;
		Vector2 directiony;
		Matrix2 covmatrix;
		CovarianceMatrix2(
			size, const_cast<Point*>(pointarray),
			&covmatrix
		);
		covmatrix.EnginVector(&directionx, &directiony);
		if (directionx.CrossProduct(directiony) > 0.0) {
			directiony.Invert();
		}
		Real xmin = 0.0;
		Real xmax = 0.0;
		Real ymin = 0.0;
		Real ymax = 0.0;
		CalculateProjectionLength(directionx, size, pointarray, &xmin, &xmax);
		CalculateProjectionLength(directiony, size, pointarray, &ymin, &ymax);
		if (ymax - ymin > xmax - xmin) {
			box_.direction_x = directiony;
			box_.direction_y = directionx;
			box_.half_width = (ymax - ymin) * 0.5;
			box_.half_height = (xmax - xmin) * 0.5;
		}
		else {
			box_.direction_x = directionx;
			box_.direction_y = directiony;
			box_.half_width = (xmax - xmin) * 0.5;
			box_.half_height = (ymax - ymin) * 0.5;
		}
		box_.center =
			directionx * ((xmin + xmax) * 0.5) +
			directiony * ((ymin + ymax) * 0.5);
	}

	void OrientedBoundingBox::CreateBounding(
		const PointSet& pointset
	) {
		PointArray pointarray(pointset);
		CreateBounding(pointarray.Size(), pointarray.Array());
	}

	void OrientedBoundingBox::CreateBounding(
		const Circle& circle
	) {
		box_.center = circle.center;
		box_.direction_x.Set(1, 0);
		box_.direction_y.Set(0, 1);
		box_.half_width = box_.half_height = circle.radius;
	}

	void OrientedBoundingBox::CreateBounding(
		const Arc& arc
	) {
		box_.direction_x = (arc.end - arc.start).NormalizeSafe();
		box_.direction_y = VectorVerticalProduct(box_.direction_x);
		box_.half_width = 0.5 * (arc.end - arc.start).DotProduct(box_.direction_x);
		box_.half_height = 0.5 * (arc.radius - (arc.start - arc.center).DotProduct(box_.direction_y));
		box_.center = arc.center + box_.direction_y * (arc.radius - box_.half_height);
	}

	void OrientedBoundingBox::CreateBounding(
		const Sector& sector
	) {
		box_.direction_x = (sector.end - sector.start).NormalizeSafe();
		box_.direction_y = VectorVerticalProduct(box_.direction_x);
		box_.half_width = 0.5 * (sector.end - sector.start).DotProduct(box_.direction_x);
		box_.half_height = 0.5 * sector.radius;
		box_.center = sector.center + box_.direction_y * box_.half_height;
	}

	void OrientedBoundingBox::CreateBounding(
		const Triangle& triangle
	) {
		CreateBounding(3, triangle.vertex);
	}

	void OrientedBoundingBox::CreateBounding(
		const Rectangle& rectangle
	) {
		box_.center = rectangle.Center();
		box_.direction_x = Vector2(1, 0);
		box_.direction_y = Vector2(0, 1);
		box_.half_width = rectangle.width * 0.5;
		box_.half_height = rectangle.height * 0.5;
	}

	void OrientedBoundingBox::CreateBounding(
		const Box& box
	) {
		box_ = box;
	}

	void OrientedBoundingBox::CreateBounding(
		const Polygon& polygon
	) {
		CreateBounding(polygon.n, polygon.vertex);
	}

	OrientedBoundingBox CombineBounding(
		const OrientedBoundingBox& obb_1,
		const OrientedBoundingBox& obb_2
	) {
		Box box1(obb_1.GetGeometry());
		Box box2(obb_2.GetGeometry());
		Point points[8] = {
			box1.center - box1.direction_x * box1.half_width -
			box1.direction_y * box1.half_height,

			box1.center + box1.direction_x * box1.half_width -
			box1.direction_y * box1.half_height,

			box1.center + box1.direction_x * box1.half_width +
			box1.direction_y * box1.half_height,

			box1.center - box1.direction_x * box1.half_width +
			box1.direction_y * box1.half_height,

			box2.center - box2.direction_x * box2.half_width -
			box2.direction_y * box2.half_height,

			box2.center + box2.direction_x * box2.half_width -
			box2.direction_y * box2.half_height,

			box2.center + box2.direction_x * box2.half_width +
			box2.direction_y * box2.half_height,

			box2.center - box2.direction_x * box2.half_width +
			box2.direction_y * box2.half_height
		};
		OrientedBoundingBox obb;
		obb.CreateBounding(8, points);
		return obb;
	}

} //namespace grid