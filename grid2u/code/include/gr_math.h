//@header	:	gr_math.h
//@abstract	:	1.Define and implement some useful math @funtion.
//@author	:	zhong
//@date		:	2018/8
//@version	:	bate0.1

//TODO(zhong) : maybe i will implement it myself.but now the point is not the math.
//				so i use the standard c++ math.

//add		:	1.add preposition declaration
//				2.add functions
//				template<typename T> T GrExpectation(int size, T* t_array);
//				template<typename T> T GrVariance(int size, T* t_array);
//@date		:	2018/8/28

#ifndef _GRID_MATH_H_
#define _GRID_MATH_H_

#include"gr_precision.h"
#include<cassert>
#include<cmath>

namespace grid
{
	//@function-->
	inline Real GrMin(Real r1, Real r2)
	{
		return r1 > r2 ? r2 : r1;
	}

	inline Real GrMax(Real r1, Real r2)
	{
		return r1 < r2 ? r2 : r1;
	}

	inline Real GrFabs(Real r) 
	{
		return std::fabs(r);
	}

	inline Real GrSqrt(Real r)
	{
		return std::sqrt(r);
	}

	inline Real GrQuickSqrtReciprocal(Real r) 
	{	//must be float
		float x = static_cast<float>(r);
		float xhalf = 0.5f * x;
		int i = *(int*)&x;
		i = 0x5f3759df - (i >> 1);
		x = *(float*)&i;
		x = x * (1.5f - xhalf * x * x);
		x = x * (1.5f - xhalf * x * x);
		x = x * (1.5f - xhalf * x * x);
		return static_cast<Real>(x);
	}

	//@para : nedd radian
	inline Real GrSin(Real r) 
	{
		return std::sin(r);
	}

	//@para : need angle
	inline Real GrSinA(Real angle)
	{
		return std::sin(AngleToRadian(angle));
	}

	//@para : nedd radian
	inline Real GrCos(Real r) 
	{
		return std::cos(r);
	}

	//@para : need angle
	inline Real GrCosA(Real angle)
	{
		return std::cos(AngleToRadian(angle));
	}

	//@para : nedd radian
	inline Real GrTan(Real r) 
	{
		return std::tan(r);
	}

	inline Real GrTanA(Real angle)
	{
		return std::tan(AngleToRadian(angle));
	}

	//@return : return radian
	inline Real GrArcSin(Real r) 
	{
		return std::asin(r);
	}

	inline Real GrArcSinA(Real r)
	{
		return RadianToAngle(std::asin(r));
	}

	//@return : return radian
	inline Real GrArcCos(Real r) 
	{
		return std::acos(r);
	}

	inline Real GrArcCosA(Real r)
	{
		return RadianToAngle(std::acos(r));
	}

	//@return : return radian
	inline Real GrArcTan(Real r) 
	{
		return std::atan(r);
	}

	inline Real GrArcTanA(Real r)
	{
		return RadianToAngle(std::atan(r));
	}

	//@abstract : mathematical expectation
	//@para	size : the given array size,must > 0
	//@para t_array :	the pointer to the given array
	//					only support type [Real] [Vector2] [Vector3]
	//@tips	: when [T = Real]
	//			the return is the average of the given array
	//		  when [T = Vector2/Vector3]
	//			the return is the average of x and y and z.
	//			respectively save in the return_vector.x or return_vector.y or return_vector.z
	template<typename T>
	T GrExpectation(int size, T* t_array)
	{
		assert(size > 0);
		assert(t_array);
		T expectation(0.0);
		for (int i = 0; i < size; ++i)
		{
			expectation += t_array[i];
		}
		return expectation / size;
	}

	//@abstract : @GrExpection
	//@tips		: only support [Real] [Vector2] [Vector3]
	template<typename T>
	T GrVariance(int size, T* t_array)
	{
		assert(size > 0);
		assert(t_array);
		T expection = GrExpectation(size, t_array);
		T variance(0.0);
		for (int i = 0; i < size; ++i)
		{
			variance += (t_array[i] - expection) * (t_array[i] - expection);
		}
		return variance / size;
	}

}//namespace grid

#endif //_GRID_MATH_H_