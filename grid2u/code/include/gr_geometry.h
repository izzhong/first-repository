//@header : gr_geometry.h
//@author : zhong
//@date	  : 2018/9

//BUG : 1.右值构造和右值赋值都有BUG
//		2.PointSet :: End() 函数有BUG 暂时不用了

#ifndef _GRID_GEOMETRY_H_
#define _GRID_GEOMETRY_H_

#include"gr_matrix.h"
#include<initializer_list>

namespace grid
{
	//@decleration
	//@abstract :	define Vector2 as 2D-Point
	//	`			os you can use a Point like a Vector2 very easily
	using Point = Vector2;
	class PointArray;
	class PointSet;

	class Geometry;
	class Circle;
	class Arc;
	class Sector;
	class Segment;
	class Ray;
	class Triangle;
	class Rectangle;
	class Box;
	class Polygon;

	//@abstract : PointArray
	//			  use it like a real array and do not care memory
	class PointArray
	{
	public:

		PointArray() {  }

		PointArray(int size)
		{
			if (size > 0)
			{
				size_ = size;
				assert(pointarray_ = new Point[size_]{});
			}
		}

		//@abstract : use a initialization list to init a point array
		//			  like this
		//			  PointArray pa{p1,p2,p3,...,pn};
		PointArray(std::initializer_list<Point> il);

		//@abstract : use another PointSet to init a point array
		PointArray(const PointSet& pointset);

		//@abstract : init a point array by copy anothor 
		PointArray(const PointArray& pointarray);

		//@abstract : init a point array by move anothor
		//PointArray(PointArray&& pointarray);

		//@abstract : assgin to a point array by copy anothor
		PointArray& operator=(const PointArray& pointarray);

		//@abstract : assgin to a point array by move anothor
		//PointArray& operator=(PointArray&& pointarray);

		//@abstract : destroy the momery
		~PointArray(){Destruct_();}

		//@abstract : destroy the point array
		void Clear() { Destruct_(); }

		//@abstract : the number of the point array
		int Size() const { return size_; }

		//@abstract : the first address of the point array
		Point* Array() const { return pointarray_; }

		//@abstract : is point array empty?
		bool Empty() const { return size_ == 0; }

		//@abstract : you can use it like a real array
		//			  like this
		//			  PointArray pa{p1,p2,...,pn};
		//			  for(int i = 0;i < pa.Size(); ++i)
		//				todo(pa[i]);
		//@return   : return the reference
		//@tips		: error index will assert
		Point& operator[](int index)
		{
			assert(index >= 0);
			assert(index < size_);
			return pointarray_[index];
		}

		//@abstract : you can use it like a real array
		//			  like this
		//			  PointArray pa{p1,p2,...,pn};
		//			  for(int i = 0;i < pa.Size(); ++i)
		//				todo(pa[i]);
		//@return   : return the value
		//@tips		: error index will assert
		Point operator[](int index) const
		{
			assert(index >= 0);
			assert(index < size_);
			return pointarray_[index];
		}

		void Transform(const Matrix3& matrix)
		{
			matrix.Transform(size_, pointarray_);
		}

	private:

		//@abstract : hold the number of the point array
		//			  default is zero
		int size_ = 0;

		//@abstract : hold the first address of the point array
		//			  default is null ptr
		Point* pointarray_ = nullptr;

	private:

		//@abstract : accocate the memory
		//			  error size will assert.
		Point* Alloc_(int size)
		{
			assert(size > 0);
			if (size == 0)
				return nullptr;
			else
				return new Point[size]{};
		}

		//@abstract : destroy the memory
		void Destruct_()
		{
			if (Empty())
				return;
			if (size_ == 1)			
				delete pointarray_;			
			else
				delete[] pointarray_;
			size_ = 0;
		}
	};

	//@abstract : using for class[PointSet]
	//@members	[point] :	the node value of [Point]
	//			[next]	:	the next node ptr
	//			[prev]	:	the previous node ptr
	struct PointNode
	{
		Point		point;
		PointNode*	next;
		PointNode*	prev;
	};

	//use point array to define fixed point set.

	//@class : PointSet
	//@abstract : unfixed point set
	class PointSet
	{
	public:
		//@abstract : define a empty point set
		PointSet() :
			size_(0), head_(nullptr)
		{	}

		//@abstract : define a point set by copy another point set
		PointSet(const PointSet& pointset);

		//@abstract : define a point set by move another point set
		//PointSet(PointSet&& pointset);

		//@abstract : define a point set by copy a given point array
		PointSet(int size, const Point* point_array);

		//@abstract	: reset the point set by copy another point set
		//			  old things will be released
		//PointSet& operator=(const PointSet& pointset);

		//@abstract	: reset the point set by move another point set
		//			  old things will be released
		PointSet& operator=(PointSet&& pointset);

		//@abstract : destruct point set and free memory
		~PointSet() { Destruct_(); }

		//@abstract : clear the point set and 
		//			  reset the point set by given point array
		PointSet& Set(int size, const Point* point_array);
		
		//@abstract : is the point set has no point
		bool Empty() const { return size_ == 0 || head_ == nullptr; }

		//@abstract : the numbers of point in this point set
		int Size() const { return size_; }

		//@abstract : release all the points
		void Clear() { Destruct_(); }

		//@abstract : the head node ptr 
		PointNode* HeadPtr() const { return head_; }

		//@abstract : the tail node ptr
		PointNode* TailPtr() const {
			if (head_)
				return head_->prev;
			else
				return nullptr;
		}

		//@abstract : only use for tranverse
		//@warning  : only use it like this
		//				for(auto work = ps.HeadPtr();ps.End(work);work = work->next or work = work->prev)
		//					todo
		bool End(const PointNode* node) const
		{
			if (Empty()) {
				return true;
			}			
			if (node == head_) {
				return (loop_count_ = !loop_count_);
			}
			else {
				return false;
			}
		}

		//@abstract : add the given point after the tail
		//@return	: the node ptr that creates newly
		PointNode* Push(const Point& point);

		//@abstract : add the given pointnode after the tail
		void Push(PointNode* pointnode);

		//@abstract : delete the head node
		void EraseHead();

		//@abstract : delte the tail node
		void EraseTail();

		//@abstract : delete the given node
		//			  if the nodeptr is null.fun do nothing.
		//			  if the nodeptr is not in the pointset,
		//			  fun will error.so use it carefully.
		void Erase(const PointNode* node);

		//@abstract : delete the node which is == the given point.
		//			  fun will only delete the first == point.
		//			  if the point set has not point == the given point.
		//			  fun will do nothing and return false,otherwise return true
		bool Erase(const Point& point);

		//@abstract : insert the given point after the given nodeptr
		//			  if the nodeptr is not in the point set.
		//			  fun will error.
		//@return	: return the newly insert nodeptr
		PointNode* InsertAfter(PointNode* at, const Point& point);

		//@abstract : insert the given point before the given nodeptr
		//			  if the nodeptr is not in the point set.
		//			  fun will error.
		//@return	: return the newly insert nodeptr
		PointNode* InsertBefore(PointNode* at, const Point& point);

		//@abstract : merge the given point set into this by copy
		PointSet& Merge(const PointSet& pointset);

		//@abstract : merge the given point set into this by move
		PointSet& Merge(PointSet&& pointset);

		//@abstract : find the nodeptr which is == given point.
		//@para [positive] : default is true,which stands for Positive sequence traversal
		PointNode* Find(const Point& point, bool positive = true);

		//@abstract : find the first nodeptr which is match the condition
		//@para [unary]	   : a callable object,which is the condition
		//		[positive] : default is true,which stands for Positive sequence traversal
		PointNode* Find(bool(*unary)(const Point&), bool positive = true);

		//@abstract : use the given callable object to call every points
		void ForEach(void(*todo)(const Point&), bool positive = true) const;

		//@abstract : use the given callable object to call every points
		void ForEach(void(*todo)(Point&), bool positive = true);

		//@abstract : use the given callable object to call every points
		void ForEach(void(*todo)(const Point*), bool positive = true) const;

		//@abstract : use the given callable object to call every points
		void ForEach(void(*todo)(Point*), bool positive = true);

		//@abstract : use the transform matrix to every points
		void Transform(const Matrix3& matrix);

	private:
		//@abstract : hold the head nodeptr
		PointNode * head_;
		//@abstract : hold the numbers of points
		int size_;
		//@abstract : support End();
		mutable bool loop_count_ = true;
	private:
		//@abstract : allocate the memory and assign
		PointNode* Alloc_(const Point& point);
		//@abstract : release all the memory
		void Destruct_();
	}; //class PointSet

	//@abstract : Geometry
	class Geometry
	{
	public:

		virtual void Transform(const Matrix3& matrix) = 0;

		virtual bool IsStandard() const = 0;
	};

	//@class : Circle
	class Circle :
		public Geometry
	{
	public:

		//@abstract : hold the center of the circle
		Point center;

		//@abstract : hold the radius of the circle
		Real radius;

	public:

		Circle() :
			center(), radius(1.0)
		{	}

		Circle(const Point& Center, Real r) :
			center(Center), radius(r)
		{	}

		Circle& Clear()
		{
			center.Clear();
			radius = 0.0;
			return *this;
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(&center);
		}

		//@abstract : define a standard circle is radius > 0
		virtual bool IsStandard() const override
		{
			return radius > 0.0;
		}

		//@abstract : the area of the circle
		Real Area() const
		{
			return MATH_CONSTANT_PI * radius * radius;
		}
	}; //class Circle

	//@class : Arc
	class Arc :
		public Geometry
	{
	public:

		//@abstract : the center of the arc
		Point center;

		//@abstract : the radius fo the arc
		Real radius;

		//@abstract : the start point at the arc
		Point start;

		//@abstract : the end point at the arc
		Point end;

	public:

		Arc() :
			center(), start(), end(), radius(0.0)
		{	}

		Arc(
			const Point& _center,
			const Point& _start,
			const Point& _end,
			Real r
		) :
			center(_center), start(_start), end(_end), radius(r)
		{	}

		Arc& Clear()
		{
			center.Clear();
			start.Clear();
			end.Clear();
			radius = 0.0;
			return *this;
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(&center);
			matrix.Transform(&start);
			matrix.Transform(&end);
		}

		//@abstract : define a standard arc
		//			  [1] radius > 0
		//			  [2] start !|| end
		//			  [3] radius = len(center,start) = len(center,end)	
		virtual bool IsStandard() const override
		{
			return radius > 0.0 &&
				IsZeroReal((start - center).LengthSquare() - radius * radius) &&
				IsZeroReal((end - center).LengthSquare() - radius * radius) &&
				!IsZeroReal(((start - center).CrossProduct(end - center)));
		}

		//@abstract : the length of the arc
		Real Length() const
		{
			Real radian = VectorRadian(start - center, end - center);
			return GrFabs(radian * radius);
		}
	}; //class Arc

	//@class : Sector
	class Sector :
		public Geometry
	{
	public:

		//@abstract : the center of the sector
		Point center;

		//@abstract : the radius fo the sector
		Real radius;

		//@abstract : the start point at the sector arc
		Point start;

		//@abstract : the end point at the sector arc
		Point end;

	public:

		Sector() :
			center(), start(), end(), radius(0.0)
		{	}

		Sector(
			const Point& _center,
			const Point& _start,
			const Point& _end,
			Real r
		) :
			center(_center), start(_start), end(_end), radius(r)
		{	}

		Sector& Clear()
		{
			center.Clear();
			start.Clear();
			end.Clear();
			radius = 0.0;
			return *this;
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(&center);
			matrix.Transform(&start);
			matrix.Transform(&end);
		}

		//@abstract : define a standard sector
		//			  [1] radius > 0
		//			  [2] start != end
		//			  [3] radius = len(center,start) = len(center,end)	
		virtual bool IsStandard() const override
		{
			return
				radius > 0.0 &&
				start != end &&
				IsZeroReal((start - center).LengthSquare() - radius * radius) &&
				IsZeroReal((end - center).LengthSquare() - radius * radius);
		}

		//@abstract : the area of the sector
		Real Area() const
		{
			Real radian = VectorRadian(start - center, end - center);
			return GrFabs(radian * radius * radius);
		}
	}; //class Sector

	   //@class : Segment
	class Segment :
		public Geometry
	{
	public:

		//@abstract : the start point of the segment
		Point start;

		//@abstract : the end point of the segment
		Point end;

	public:

		Segment() :
			start(), end(1.0, 0.0)
		{	}

		Segment(const Point& start_point, const Point& end_point) :
			start(start_point), end(end_point)
		{	}

		Segment& Clear()
		{
			start.Clear();
			end.Clear();
			return *this;
		}

		Segment& Set(const Point& start_, const Point& end_)
		{
			start = start_;
			end = end_;
			return *this;
		}

		Segment& Set(const Segment& segment)
		{
			return *this = segment;
		}

		//@abstract : define a standard segment is two different points
		virtual bool IsStandard() const override
		{
			return start != end;
		}

		//@abstract : the length of the segment
		Real Length() const
		{
			return (end - start).Length();
		}

		//@abstract : the square of the segment length
		Real LengthSquare() const
		{
			return (end - start).LengthSquare();
		}

		//@abstract : the normalize normal vector of the segment
		//			  if the segment is not standard,the return is zero vector.
		//@para [clock_wise] : true is clock_wise;false is unti clock wise
		Vector2 Normal(bool clock_wise = true) const;

		//@abstract : the normalize direction vector of the segment by start to end
		//			  if the segment is not standard,the return is zero vector.
		Vector2 Direction() const
		{
			if (!IsStandard()) {
				return Vector2(0.0, 0.0);
			}
			return (end - start).Normalize();
		}

		//@abstract : the center point the segment
		Point Center() const
		{
			return (start + end) * 0.5;
		}

		//@abstract : transform the start and end point
		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(&start);
			matrix.Transform(&end);
		}

	}; //class Segment

	   //@class : Ray
	class Ray :
		public Geometry
	{
	public:

		//@abstract : the origin of the ray
		Point origin;

		//@abstract : the direction of the ray
		Vector2 direction;

	public:

		Ray() :
			origin(), direction(1.0, 0.0)
		{	}

		Ray(const Point& _origin, const Vector2& _direction) :
			origin(_origin), direction(_direction)
		{	}

		Ray& Clear()
		{
			origin.Clear();
			direction.Clear();
			return *this;
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			Point P(origin + direction);
			matrix.Transform(&origin);
			matrix.Transform(&P);
			direction = (P - origin);
		}

		//@abstract : a standard ray's direction vector is normalization vector.
		virtual bool IsStandard() const	override {
			return direction.IsNormalization();
		}
	}; //class Ray

	//@class : Line
	class Line :
		public Geometry
	{
	public:

		Point origin;
		Vector2 normal;

	public:

		Line(const Point& origin_,const Vector2& normal_) :
			origin(origin_),normal(normal_)
		{	}

		Line(const Segment& segment) :
			origin(segment.start),
			normal(VectorVerticalProduct((segment.end - segment.start)).NormalizeSafe())
		{	}

		Line(const Ray& ray) :
			origin(ray.origin),
			normal(VectorVerticalProduct(ray.direction))
		{	}

		bool IntersectionPoint(const Line& line,Point* point) const
		{
			assert(point);
			//对于两条直线的交点可以使用矩阵乘法
			Matrix2 m{
				normal.x,normal.y,
				line.normal.x,line.normal.y
			};
			if (m.IsInvertible()) {
				point->Set(
					normal.DotProduct(origin), 
					line.normal.DotProduct(line.origin)
				);
				point->Set(m.Invert()*(*point));
				return true;
			}
			else{
				return false;
			}			
		}

		Point Intersection(const Line& line) const 
		{
			Point point;
			Matrix2 m{
				normal.x,normal.y,
				line.normal.x,line.normal.y
			};
			assert(m.IsInvertible());
			point.Set(
				normal.DotProduct(origin),
				line.normal.DotProduct(line.origin)
			);
			point.Set(m.Invert()*(point));
			return point;
		}

		bool IsIntersectant(const Line& line) const
		{
			return !IsVectorParallel(normal, line.normal);
		}

		//Ax + By = C
		Vector3 LinearEquation() const
		{
			return Vector3(normal.x, normal.y, 
				origin.x * normal.x + origin.y * normal.y);
		}

		virtual void Transform(const Matrix3& matrix) override
		{

		}

		virtual bool IsStandard() const override
		{
			return true;
		}

	}; //class Line

	

	

	//@class : Triangle
	class Triangle :
		public Geometry
	{
	public:

		//@abstract : hold the three vertexs of a triangle
		Point vertex[3];

	public:

		Triangle() :
			vertex{  }
		{	}

		Triangle(
			const Point& vertex_1,
			const Point& vertex_2,
			const Point& vertex_3
		) :
			vertex{vertex_1,vertex_2,vertex_3}
		{	}

		Triangle& Clear()
		{
			vertex[0].Clear();
			vertex[1].Clear();
			vertex[2].Clear();
			return *this;
		}

		Triangle& Set(
			const Point& vertex_0,
			const Point& vertex_1,
			const Point& vertex_2
		) {
			vertex[0] = vertex_0;
			vertex[1] = vertex_1;
			vertex[2] = vertex_2;
			return *this;
		}

		bool ClockWise() const
		{
			return (vertex[1] - vertex[0]).CrossProduct(vertex[2] - vertex[1]) > 0.0;
		}

		//@abstract : define a standard triangle is non-collinear three points
		virtual bool IsStandard() const override
		{
			return !IsZeroReal(
				(vertex[1] - vertex[0]).CrossProduct(vertex[2] - vertex[0])
			);
		}

		//@abstract : get the side segment
		//@para [index] : range[0,1,2].
		//				  side[0] is segment(vertex[0],vertex[1])
		//				  side[1] is segment(vertex[1],vertex[2])
		//				  side[2] is segment(vertex[2],vertex[0])
		//				  other input will assert
		Segment Side(int index) const
		{
			assert(index >= 0 && index < 3);
			return Segment(vertex[index], vertex[(index + 1) % 3]);
			//[index + 1 > 3 ? 0 : index + 1] which is faster?
		}

		//@abstract : get the normalization side normal vector
		//@para [index] : range[0,1,2].
		//				  side normal[0] is segment(vertex[0],vertex[1]).normal()
		//				  side normal[1] is segment(vertex[1],vertex[2]).normal()
		//				  side normal[2] is segment(vertex[2],vertex[0]).normal()
		//				  other input will assert
		//				  if the side segment is not standard,result is zero vector.
		Vector2 SideNormal(int index) const;

		//@abstract : the center point of the triangle
		Point Center() const
		{
			return (vertex[1] + vertex[2] + vertex[0]) * (0.333333333333333333333333);
		}

		//@abstract : the area of the triangle
		Real Area() const
		{
			return 0.5 * GrFabs((vertex[1] - vertex[0]).CrossProduct(vertex[2] - vertex[0]));
		}

		//@abstract : transform the triangle by the matrix
		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(3, vertex);
		}
	}; //class Triangle

	//@class : Rectangle
	class Rectangle :
		public Geometry
	{
	public:

		//@abstract : hold the left x
		Real x;

		//@abstract : hold the top y
		Real y;

		//@abstract : hold the width of the rectangle
		Real width;

		//@abstract : hold the height of the rectangle
		Real height;

	public:

		Rectangle() :
			x(0.0), y(0.0), width(0.0), height(0.0)
		{	}

		Rectangle(Real _x, Real _y, Real _width, Real _height) :
			x(_x), y(_y), width(_width), height(_height)
		{	}

		//@abstract : use two main dialog point to create a rectangle 
		Rectangle(const Point& lt_point, const Point& rb_point) :
			x(lt_point.x), y(lt_point.y),
			width(rb_point.x - lt_point.x),
			height(rb_point.y - lt_point.y)
		{
			if (width < 0 || height < 0) {
				x = rb_point.x;
				y = rb_point.y;
				width = -width;
				height = -height;
			}
		}

		//@abstract : use leading diagonal to define a rectangle
		Rectangle(const Segment& segment) :
			x(std::fmin(segment.start.x, segment.end.x)),
			y(std::fmin(segment.start.y, segment.end.y)),
			width(std::fabs(segment.start.x - segment.end.x)),
			height(std::fabs(segment.start.y - segment.end.y))
		{	}


		Box ToBox() const;

		Rectangle& Clear()
		{
			x = y = width = height = 0.0;
			return *this;
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			//用主对角线模拟可能更好一点??
			Point xy(x, y);
			Point wh(width, height);
			matrix.Transform(&xy);
			matrix.Transform(&wh);
			x = xy.x;
			y = xy.y;
			width = wh.x;
			height = wh.y;
		}

		//@abstract : define a standard rectangle is a non-zero rectangle
		virtual bool IsStandard() const override
		{
			return !(IsZeroReal(width) || IsZeroReal(height));
		}

		//@abstract : the center point of the rectangle
		Point Center() const
		{
			return Point(x + 0.5 * width, y + 0.5 * height);
		}

		//@abstract : the area of the rectangle
		Real Area() const
		{
			return width * height;
		}
	}; //class Rectangle

	//@class : OOB
	class Box :
		public Geometry
	{
	public:

		//@abstract : the center point of the box
		Point center;

		//@abstract : the direction of width
		Vector2 direction_x; 

		//@abstract : the direction of height
		Vector2 direction_y; 

		//@abstract : the half width in direction x
		Real half_width;

		//@abstract : the half height in direction y
		Real half_height;
		
	public:

		Box() :
			center(), 
			direction_x(1.0, 0.0), 
			direction_y(0.0, 1.0),
			half_width(1.0), 
			half_height(1.0)
		{	}

		Box(const Point& _center,
			const Vector2& dir_x,
			const Vector2& dir_y,
			Real h_width,
			Real h_height
		) :
			center(_center),
			direction_x(dir_x),
			direction_y(dir_y),
			half_width(h_width),
			half_height(h_height)
		{	}

		Box& Clear()
		{
			center.Clear();
			direction_x.Clear();
			direction_y.Clear();
			half_width = half_height = 0.0;
			return *this;
		}

		//@abstract : define a standard box
		//			  [1] direction_x and direction_y is normalization
		//			  [2] area is not zero
		virtual bool IsStandard() const override
		{
			return 
				direction_x.IsNormalization() && 
				direction_y.IsNormalization() && 
				half_height != 0.0 && 
				half_width != 0.0;
		}

		//@abstract : the area of the box
		Real Area() const
		{
			return 4.0 * half_width * half_height;
		}

		//@abstract : get the rectangle by rotate the box at center
		Rectangle ToRectangle() const
		{
			return Rectangle(
				center.x - half_width,
				center.y - half_height,
				2.0 * half_width,
				2.0 * half_height
			);
		}

		//@abstract : get the rectangle by the rotate box at center
		//			  and the rotate matrix
		void ToRectangle(Rectangle* rectangle, Matrix3* rotate_matrix) const;

		//@abstract : transform a box
		//			  when scale and shear may be failed
		virtual void Transform(const Matrix3& matrix)override;
	}; //class Box

	//@class : Polygon
	class Polygon :
		public Geometry
	{
	public:

		//@abstract : hold the fixed vertex point array
		const int n;
		Point * vertex;

	public:

		//@abstract : new when create
		Polygon(int size) :
			n(size),
			vertex(nullptr)		
		{
			assert(n > 0);
			assert(vertex = new Point[n]{});
		}

		//@abstract : new when create
		Polygon(int size, Point* point_array) :	
			n(size),
			vertex(nullptr)
		{
			assert(n > 0);
			assert(vertex = new Point[n]);
			for (int i = 0; i < size; ++i)
				vertex[i] = point_array[i];
		}

		Polygon(const PointSet& pointset) :
			n(pointset.Size()),
			vertex(nullptr)
		{
			if (!pointset.Empty()) {
				assert(vertex = new Point[n]{});
				PointNode* work = pointset.HeadPtr();
				int index = 0;
				do {
					vertex[index++] = work->point;
					work = work->next;
				} while (work != pointset.HeadPtr());
			}
		}

		Polygon(std::initializer_list<Point> il) :
			n(il.size()),
			vertex(nullptr)			
		{
			assert(vertex = new Point[n]);
			int i = 0;
			for (auto it = il.begin(); it != il.end(); ++it)
			{
				vertex[i++] = *it;
			}
		}

		Polygon(const Circle& circle, int precision = 0) :
			n(4 * precision + 12),
			vertex(nullptr)			
		{	//优化方法 计算出圆形的1/4之后 其余的点可以通过镜像生成
			assert(n >= 12);
			assert(vertex = new Point[n]{});
			Real delta_angle = 360.0 / n;
			Matrix3 rotate;
			RotateMatrix(Vector2(0, 0), delta_angle, &rotate);
			Vector2 radius_vector(1, 0);
			vertex[0] = circle.center + radius_vector * circle.radius;
			for (int i = 1; i < n; ++i)
			{	
				rotate.Transform(&radius_vector);
				vertex[i] = circle.center + radius_vector * circle.radius;
			}
		}

		Polygon(const Arc& arc, int precision = 0) :
			n(4 * precision + 12),
			vertex(nullptr)
		{
			assert(n >= 2);
			assert(vertex = new Point[n]{});
			Real delta_angle = VectorAngle(arc.start, arc.end) / (n - 1);
			Matrix3 rotate;
			RotateMatrix(Vector2(0, 0), delta_angle, &rotate);
			Vector2 radius_vector(arc.start - arc.center);
			vertex[0] = arc.center + radius_vector;
			for (int i = 1; i < n; ++i)
			{
				rotate.Transform(&radius_vector);
				vertex[i] = arc.center + radius_vector;
			}
		}

		Polygon(const Sector& sector, int precision = 0) :
			n(1 + 4 * precision + 12),
			vertex(nullptr)
		{
			assert(n >= 3);
			assert(vertex = new Point[n]{});
			vertex[0] = sector.center;
			Real delta_angle = VectorAngle(sector.start, sector.end) / (n - 2);
			Matrix3 rotate;
			RotateMatrix(Vector2(0, 0), delta_angle, &rotate);
			Vector2 radius_vector(sector.start - sector.center);
			vertex[1] = sector.center + radius_vector;
			for (int i = 2; i < n; ++i)
			{
				rotate.Transform(&radius_vector);
				vertex[i] = sector.center + radius_vector;
			}
		}

		Polygon(const Segment& segment) :
			n(2),
			vertex(nullptr)
		{
			assert(vertex = new Point[n]{});
			vertex[0] = segment.start;
			vertex[1] = segment.end;
		}

		Polygon(const Ray& ray, Real length = 1.0):
			n(2),
			vertex(nullptr)
		{
			assert(vertex = new Point[n]{});
			vertex[0] = ray.origin;
			vertex[1] = ray.origin + ray.direction * length;
		}

		Polygon(const Triangle& triangle) :
			n(3),
			vertex(nullptr)
		{
			assert(vertex = new Point[n]{});
			vertex[0] = triangle.vertex[0];
			vertex[1] = triangle.vertex[1];
			vertex[2] = triangle.vertex[2];
		}

		Polygon(const Rectangle& rectangle) :
			n(4),
			vertex(nullptr)
		{
			assert(vertex = new Point[n]{});
			vertex[0].Set(rectangle.x, rectangle.y);
			vertex[1].Set(rectangle.x + rectangle.width, rectangle.y);
			vertex[2].Set(rectangle.x + rectangle.width, rectangle.y + rectangle.height);
			vertex[3].Set(rectangle.x, rectangle.y + rectangle.height);
		}

		Polygon(const Box& box) :
			n(4),
			vertex(nullptr)
		{
			assert(vertex = new Point[n]{});
			vertex[0] = 
				box.center -
				box.direction_x * box.half_width - 
				box.direction_y * box.half_height;
			vertex[1] =
				box.center +
				box.direction_x * box.half_width -
				box.direction_y * box.half_height;
			vertex[2] =
				box.center +
				box.direction_x * box.half_width +
				box.direction_y * box.half_height;
			vertex[3] =
				box.center -
				box.direction_x * box.half_width +
				box.direction_y * box.half_height;
		}

		Polygon& Clear()
		{
			for (int i = 0; i < n; ++i)
			{
				vertex[i].Clear();
			}
			return *this;
		}

		Polygon& operator<<=(const Point& point)
		{
			index_ = 0;
			vertex[index_] = point;
			return *this;
		}

		Polygon& operator,(const Point& point)
		{
			assert(++index_ < n);
			vertex[index_] = point;
			return *this;
		}

		~Polygon()
		{
			if (n > 1)
				delete[] vertex;
			else
				delete vertex;
		}

		//@abstract : define a standard polygon is convex polygon
		virtual bool IsStandard() const override;

		//@abstract : the center point of the polygon
		Point Center() const;

		//@abstract : the area of the polygon
		Real Area() const;

		//@abstract : the side segment of the polygon
		//@para [index] : range [0,n-1]
		//				  side[0] is segment(vertex[0],vertex[1])
		//					...
		//				  side[n-1] is segment(vertex[n-1],vertex[0])
		//				  other input will assert
		//				  unstandard geometry will return zero vector
		Segment Side(int index) const
		{
			assert(index >= 0 && index < n);
			return Segment(vertex[index], vertex[(index + 1) % n]);
		}

		//@abstract : the side segment of the polygon
		//@para [index] : range [0,n-1]
		//				  side normal[0] is segment(vertex[0],vertex[1]),normal
		//					...
		//				  side normal[n-1] is segment(vertex[n-1],vertex[0]),normal
		//				  other input will assert
		//				  unstandard geometry will return zero vector
		Vector2 SideNormal(int index) const
		{
			assert(index >= 0 && index < n);
			Triangle triangle(
				vertex[index], vertex[(index + 1) % n], vertex[(index + 2) % n]
			);
			return triangle.SideNormal(0);
		}

		//@abstract : transform the polygon
		virtual void Transform(const Matrix3& matrix) override
		{
			matrix.Transform(n, vertex);
		}
	private:
		//@abstract : support operator<<= and operator,
		int index_ = 0;
	}; //class Polygon

	struct PolygonNode
	{
		Polygon polygon;
		PolygonNode* next;
		PolygonNode* prev;
	};

	//@class : Combined Polygon
	class CombinedPolygon :
		public Geometry
	{
	public:
	};

} //namespace grid

#endif //_GRID_GEOMETRY_H_