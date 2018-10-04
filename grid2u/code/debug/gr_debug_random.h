//@header : gr_debug_random.h
//@author : zhong
//@abstract : 用于随机生成测试元素

#ifndef _GRID_DEBUG_RANDOM_H_
#define _GRID_DEBUG_RANDOM_H_

#include<gr_geometry.h>
#include<gr_precision.h>
#include<ctime>
#include<random>

namespace grid
{
	namespace grid_debug
	{
		void RandomRealArray(
			int size, Real* real_array,
			Real min = 0.0, Real max = 1.0
		) {
			std::default_random_engine e(static_cast<unsigned int>(time(0)));
			std::uniform_real_distribution<Real> u(min, max);
			for (int i = 0; i < size; ++i)
			{
				real_array[i] = u(e);
			}
		}

		Vector2 RandomVector2(
			Real x_min = 0.0, Real x_max = 1.0,
			Real y_min = 0.0, Real y_max = 1.0
		) {
			std::default_random_engine e(static_cast<unsigned int>(time(0)));
			std::uniform_real_distribution<Real> ux(x_min, x_max);
			std::uniform_real_distribution<Real> uy(y_min, y_max);
			return Vector2(ux(e), uy(e));
		}

		void RandomVector2Array(
			int size, Vector2* vector_array,
			Real x_min = 0.0, Real x_max = 1.0,
			Real y_min = 0.0, Real y_max = 1.0
		) {
			std::default_random_engine e(static_cast<unsigned int>(time(0)));
			std::uniform_real_distribution<Real> ux(x_min, x_max);
			std::uniform_real_distribution<Real> uy(y_min, y_max);
			for (int i = 0; i < size; ++i)
			{
				vector_array[i].x = ux(e);
				vector_array[i].y = uy(e);
			}
		}

		//@class RandomGeometry
		class RandomGeometry
		{
		public:

			RandomGeometry(
				unsigned seed,
				Real x_min = 0.0, Real x_max = 1.0,
				Real y_min = 0.0, Real y_max = 1.0
			) :
				e_(seed),
				urx_(x_min,x_max),
				ury_(y_min,y_max)		
			{	}

			RandomGeometry(
				Real x_min = 0.0, Real x_max = 1.0,
				Real y_min = 0.0, Real y_max = 1.0
			) :
				e_(static_cast<unsigned int>(time(0))),
				urx_(x_min, x_max),
				ury_(y_min, y_max)
			{	}

			Point RandomPoint_() {
				return Point(urx_(e_), ury_(e_));
			}

			Segment RandomSegment() {
				return Segment(RandomPoint_(), RandomPoint_());
			}

			Ray RandomRay() {
				return Ray(RandomPoint_(), RandomPoint_().NormalizeSafe());
			}

			Line RandomLine() {
				return Line(RandomPoint_(), RandomPoint_().NormalizeSafe());
			}

			Triangle RandomTriangle() {
				return Triangle(RandomPoint_(), RandomPoint_(), RandomPoint_());
			}

			Rectangle RandomRectangle() {
				return Rectangle(RandomSegment());
			}

			Box RandomBox() {
				Vector2 dir(RandomPoint_().NormalizeSafe());
				return Box(RandomPoint_(), dir, VectorVerticalProduct(dir), urx_(e_), ury_(e_));
			}

			Arc RandomArc() {		
				Point center(RandomPoint_());
				Point start(RandomPoint_());
				Real radius = (start - center).Length();
				Point end(RandomPoint_());
				end = center + end.NormalizeSafe() * radius;
				return Arc(center, start, end, radius);
			}

			Sector RandomSector() {
				Point center(RandomPoint_());
				Point start(RandomPoint_());
				Real radius = (start - center).Length();
				Point end(RandomPoint_());
				end = center + end.NormalizeSafe() * radius;
				return Sector(center, start, end, radius);
			}

			Circle RandomCircle() {
				Point center(RandomPoint_());
				Point point(RandomPoint_());
				Real radius = (point - center).Length();
				return Circle(center, radius);
			}

		private:

			std::default_random_engine e_;
			std::uniform_real_distribution<Real> urx_;
			std::uniform_real_distribution<Real> ury_;

		}; // class RandomGeometry

	} //namespace grid_debug

} //namespace grid

#endif //_GRID_DEBUG_RANDOM_H_
