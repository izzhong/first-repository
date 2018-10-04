//@header	:	gr_vector.h
//@abstract	:	1.Define @class Vector2 and @class Vector3.
//				2.Implement some math @function relate to Vector2 and Vector3.
//@author	:	zhong
//@date		:	2018/8
//@version	:	bate1.0

//@add		:	add operator overloading in [Vector2]&[Vector3]
//				Vector2 operator*(const Vector2& right) const;
//				Vector2& operator*=(const Vector2& right);
//				Vector3 operator*(const Vector3& right) const
//				Vector3& operator*=(const Vector3& right)
//				Real Covariance(int size, Vector2* vector_array);
//				Real CorrelationCoefficient(int size, Vector2* vector_array);
//				Vector3 Covariance(int size, Vector3* vector_array)
//@date		:   2018/8/28

#ifndef _GRID_VECTOR_H_
#define _GRID_VECTOR_H_

#include"gr_math.h"

namespace grid
{	
	//@class-->
	
	//@abstract	: use it like 2d-vector : (x,y).
	class Vector2
	{
	public:	
		//the x of the 2d-vector.
		Real x;
		//the y of the 2d-vector.
		Real y;
	public:
		//default 2d-vector is (0.0,0.0).
		Vector2() :
			x(0.0), y(0.0)
		{	}

		//default y is 0.0
		explicit Vector2(Real _x) :
			x(_x),y(0.0)
		{	}

		Vector2(Real _x, Real _y) :
			x(_x), y(_y)
		{	}

		//set 2d-vector's value
		Vector2& Set(Real _x, Real _y) 
		{
			x = _x;
			y = _y;
			return *this;
		}

		//use another 2d-vector to set this
		Vector2& Set(const Vector2& vector) 
		{
			x = vector.x;
			y = vector.y;
			return *this;
		}

		//@abstract : clear a 2d vector to (0,0)
		Vector2& Clear()
		{
			x = 0.0;
			y = 0.0;
			return *this;
		}

		//invert a 2d-vector
		Vector2& Invert()
		{
			x = -x;
			y = -y;
			return *this;
		}

		//judge a 2dvector is a zero vecotr ro not
		bool IsZero() const
		{	//remove precision error
			return IsZeroReal(x) && IsZeroReal(y);
		}

		//calculate the length
		Real Length() const
		{
			return GrSqrt(x * x + y * y);
		}

		//calculate the square of the length
		Real LengthSquare() const
		{
			return x * x + y * y;
		}

		//calculate the reciprocal of the length
		Real LengthReciprocal() const
		{	//zero vector can not calculate this.
			assert(!IsZero());
			//return GrQuickSqrtReciprocal(x * x + y * y);
			return 1.0 / GrSqrt(x * x + y * y);
		}
	
		//2d-vector dot product
		Real DotProduct(const Vector2& right) const
		{
			return x * right.x + y * right.y;
		}

		//2d-vector corss product
		Real CrossProduct(const Vector2& right) const
		{
			return x * right.y - y * right.x;
		}

		//normalize a 2d-vector to let its length change to 1.
		//if you normalize a zero vector,function will assert.
		Vector2& Normalize()
		{
			Real length_reciprocal = LengthReciprocal();
			x *= length_reciprocal;
			y *= length_reciprocal;
			return *this;
		}

		Vector2& NormalizeSafe()
		{
			if (!IsZero())
				Normalize();
			return *this;
		}

		bool IsNormalization() const
		{
			return IsZeroReal(x * x + y * y - 1.0);
		}

		//@operator overloading-->>

		Vector2 operator+(const Vector2& right) const
		{	//2d-vector plus : *this + right
			return Vector2(x + right.x, y + right.y);
		}

		Vector2 operator-(const Vector2& right) const
		{	//vector subtract : *this - right
			return Vector2(x - right.x, y - right.y);
		}

		Vector2 operator-() const
		{	
			return Vector2(-x, -y);
		}

		Vector2 operator*(Real r) const
		{	//vector scalar multi : *this * r
			return Vector2(x * r, y * r);
		}

		//@abstract : component multiplication
		Vector2 operator*(const Vector2& right) const
		{
			return Vector2(x * right.x, y * right.y);
		}

		Vector2 operator/(Real r) const
		{	//devide by 0.0 is error.which will calls an assert.
			assert(r != 0.0);
			return Vector2(x / r, y / r);
		}

		Vector2& operator+=(const Vector2& right)
		{	//vector plus_equal : *this = *this + right
			x += right.x;
			y += right.y;
			return *this;
		}

		Vector2& operator-=(const Vector2& right)
		{	//vector subtract_equal : *this = *this - right
			x -= right.x;
			y -= right.y;
			return *this;
		}

		Vector2& operator*=(Real r)
		{	//vector scalar_multi_equal : *this = *this * r
			x *= r;
			y *= r;
			return *this;
		}

		//@abstract : compoment multiplication
		Vector2& operator*=(const Vector2& right)
		{
			x *= right.x;
			y *= right.y;
			return *this;
		}

		Vector2& operator/=(Real r)
		{	//devide by 0.0 is error.which will calls an assert.
			assert(r != 0.0);
			x /= r;
			y /= r;
			return *this;
		}

		bool operator==(const Vector2& right) const
		{	//remove the precision error
			return 
				IsZeroReal(x - right.x) && 
				IsZeroReal(y - right.y);
		}

		bool operator!=(const Vector2& right) const
		{
			return !(*this == right);
		}
	};//class Vector2

	//@abstract	: use it like 3d-vector : (x,y,z)
	class Vector3
	{
	public:
		//the x of the 3d-vector
		Real x;
		//the y of the 3d-vector
		Real y;
		//the z of the 3d-vector
		Real z;
	public:
		//default 3d-vector is (0.0,0.0,0.0)
		Vector3() :
			x(0.0), y(0.0), z(0.0) 
		{	}

		explicit Vector3(Real _x) :
			x(_x), y(0.0), z(0.0)
		{	}

		explicit Vector3(Real _x,Real _y) :
			x(_x),y(_y),z(0.0)
		{	}

		Vector3(Real _x, Real _y, Real _z) :
			x(_x), y(_y), z(_z) 
		{	}

		//set 3d-vector's value
		Vector3& Set(Real _x, Real _y,Real _z)
		{
			x = _x;
			y = _y;
			z = _z;
			return *this;
		}

		Vector3 Set(const Vector3 vector)
		{
			x = vector.x;
			y = vector.y;
			z = vector.z;
			return *this;
		}

		//@abstract : clear a 3d vector to (0,0,0)
		Vector3& Clear()
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;
			return *this;
		}

		//invert a 3d-vector
		Vector3& Invert()
		{
			x = -x;
			y = -y;
			z = -z;
			return *this;
		}

		//judge a 3dvector is a zero vecotr ro not
		bool IsZero() const
		{	//remove precision error
			return IsZeroReal(x) && IsZeroReal(y) && IsZeroReal(z);
		}

		//calculate the length
		Real Length() const
		{
			return GrSqrt(x * x + y * y + z * z);
		}

		//calculate the square of the length
		Real LengthSquare() const
		{
			return x * x + y * y + z * z;
		}

		//calculate the reciprocal of the length
		Real LengthReciprocal() const
		{	//zero vector can not calculate this.
			assert(!IsZero());
			return GrQuickSqrtReciprocal(x * x + y * y + z * z);
		}

		//3d-vector dot product
		Real DotProduct(const Vector3& right) const
		{
			return x * right.x + y * right.y + z * right.z;
		}

		//3d-vector corss product
		Vector3 CrossProduct(const Vector3& right) const
		{
			return Vector3(
				y * right.z - z * right.y,
				z * right.x - x * right.z,
				x * right.y - y * right.x
			);
		}

		//normalize a 3d-vector to let its length change to 1.
		//if you normalize a zero vector,function will assert.
		Vector3& Normalize()
		{
			Real length_reciprocal = LengthReciprocal();
			x *= length_reciprocal;
			y *= length_reciprocal;
			z *= length_reciprocal;
			return *this;
		}

		Vector3& NormalizeSafe()
		{
			if (!IsZero())
				Normalize();
			return *this;
		}

		bool IsNormalization() const
		{
			return IsZeroReal(x * x + y * y + z * z - 1.0);
		}

		//@operator overloading-->>

		Vector3 operator+(const Vector3& right) const
		{	//3d-vector plus : *this + right
			return Vector3(x + right.x, y + right.y, z + right.z);
		}

		Vector3 operator-(const Vector3& right) const
		{	//vector subtract : *this - right
			return Vector3(x - right.x, y - right.y, z - right.z);
		}

		Vector3 operator-() const
		{
			return Vector3(-x, -y, -z);
		}

		Vector3 operator*(Real r) const
		{	//vector scalar multi : *this * r
			return Vector3(x * r, y * r, z*r);
		}

		//@abstract : component multiplication
		Vector3 operator*(const Vector3& right) const
		{
			return Vector3(
				x * right.x, y * right.y, z * right.z
			);
		}

		Vector3 operator/(Real r) const
		{	//devide by 0.0 is error.which will calls an assert.
			assert(r != 0.0);
			return Vector3(x / r, y / r, z / r);
		}

		Vector3& operator+=(const Vector3& right)
		{	//vector plus_equal : *this = *this + right
			x += right.x;
			y += right.y;
			z += right.z;
			return *this;
		}

		Vector3& operator-=(const Vector3& right)
		{	//vector subtract_equal : *this = *this - right
			x -= right.x;
			y -= right.y;
			z -= right.z;
			return *this;
		}

		Vector3& operator*=(Real r)
		{	//vector scalar_multi_equal : *this = *this * r
			x *= r;
			y *= r;
			z *= r;
			return *this;
		}

		//@abstract : component multiplication
		Vector3& operator*=(const Vector3& right)
		{
			x *= right.x;
			y *= right.y;
			z *= right.z;
			return *this;
		}

		Vector3& operator/=(Real r)
		{	//devide by 0.0 is error.which will calls an assert.
			assert(r != 0.0);
			x /= r;
			y /= r;
			z /= r;
			return *this;
		}

		bool operator==(const Vector3& right) const
		{	//remove the precision error
			return
				IsZeroReal(x - right.x) &&
				IsZeroReal(y - right.y) &&
				IsZeroReal(z - right.z);
		}

		bool operator!=(const Vector3& right) const
		{
			return !(*this == right);
		}
	};

	//@function-->

	//@abstract : use to get a 3d-vector length
	//as for 2d-vector and 3d-vector 
	//the CrossProduct's return type is different
	//so write the template is difficult
	//use this function can make this easiler
	inline Real _real(const Vector3& vector)
	{//
		return _real(
			vector.z < 0.0 ? -vector.Length() : vector.Length()
		);
	}

	//@abstract : calculate something emm...
	//inline Real VectorTripleScalar(
	//	const Vector2& A,
	//	const Vector2& B,
	//	const Vector2& P
	//) {
	//	Vector2 normal(B - A);
	//	return 0.0;
	//}

	//@abstract : check the two vector is parallel or not
	//@para		: you must use the Vector2 or Vector3
	template<typename Vector>
	inline bool IsVectorParallel(const Vector& left, const Vector& right)
	{
		return _real(left.CrossProduct(right)) == 0.0;
	}

	//@abstract : check the two vector is vertical or not
	//@para		: you must use the Vector2 or Vector3
	template<typename Vector>
	inline bool IsVectorVertical(
		const Vector& left, const Vector& right
	) {
		return IsZeroReal(left.DotProduct(right),1e-5);
	}

	//@abstract : get a arbitrarily vertical normal vector
	//@para		: you must use the Vector2 or Vector3
	//@tip		: not normalize
	template<typename Vector>
	inline Vector VectorVerticalProduct(const Vector& vector)
	{
		return Vector(-vector.y, vector.x);
	}

	//@abstract : vector triple procuct : a x b x c
	//@para		: you must use the Vector2 or Vector3
	template<typename Vector>
	inline Vector VectorTripleProduct(
		const Vector& a, 
		const Vector& b, 
		const Vector& c
	){	//(a x b) x c = -a*(c.b) + b*(c.a)
		return -a * c.DotProduct(b) + b * c.DotProduct(a);
	}
	
	//@abstract : calculate the radian of the start-vector to end-vector : [0,pi]
	Real VectorRadian(const Vector2& start, const Vector2& end);

	//@abstract : calculate the angle of the start to end
	Real VectorAngle(const Vector2& start, const Vector2& end);

	//@abstract : 
	Real Covariance(int size, Vector2* vector_array);

	//@abstract : system API
	//			  x = cov(x,y)
	//			  y = cov(y,z)
	//			  z = cov(z,x)
	Vector3 Covariance(int size, Vector3* vector_array);

	//@abstract : 
	Real CorrelationCoefficient(int size, Vector2* vector_array);

}//namespace grid

#endif //_GRID_VECTOR_H_