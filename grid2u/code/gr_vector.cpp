//@source
//@date		:	2018/8

#include"include/gr_vector.h"

namespace grid
{
	Real Covariance(int size, Vector2* vector_array)
	{
		assert(size > 0);
		assert(vector_array);
		Vector2 expection = GrExpectation(size, vector_array);
		Real covariance = 0.0;
		for (int i = 0; i < size; ++i)
		{
			covariance +=
				(vector_array[i].x - expection.x) *
				(vector_array[i].y - expection.y);
		}
		return covariance / static_cast<Real>(size);
	}

	Vector3 Covariance(int size, Vector3* vector_array)
	{
		assert(size > 0);
		assert(vector_array);
		Vector3 expection = GrExpectation(size, vector_array);
		Vector3 covariance;
		for (int i = 0; i < size; ++i)
		{
			covariance.x +=
				(vector_array[i].x - expection.x) *
				(vector_array[i].y - expection.y);
			covariance.y +=
				(vector_array[i].y - expection.y) *
				(vector_array[i].z - expection.z);
			covariance.z +=
				(vector_array[i].z - expection.z) *
				(vector_array[i].x - expection.x);
		}
		return covariance /= static_cast<Real>(size);
	}

	Real CorrelationCoefficient(int size, Vector2* vector_array)
	{
		assert(size > 0);
		assert(vector_array);
		Real covariance = Covariance(size, vector_array);
		Vector2 variance = GrVariance(size, vector_array);
		Real r = GrSqrt(variance.x) * GrSqrt(variance.y);
		//error : devide by zero 
		if (IsZeroReal(r)) {
			return 0.0;
		}
		else {
			return covariance / r;
		}
	}

	//函数计算失误的原因大部分都是浮点数非常小导致的
	//我们可以检测浮点数的模长 如果非常小的话 可以进行扩大
	//不一定有用 最后优化的时候再做测试
	Real VectorRadian(const Vector2& start, const Vector2& end)
	{
		if (start.IsZero() || end.IsZero())
			return 0.0;
		Vector2 a(start);
		Vector2 b(end);
		if (a.LengthSquare() < 1.0) {
			a *= 10000.0;
		}
		if (b.LengthSquare() < 1.0) {
			b *= 10000.0;
		}
		if (IsVectorParallel(a, b)) {
			if (a.DotProduct(b) > 0.0) {
				return 0.0;
			}
			else {
				return MATH_CONSTANT_PI;
			}
		}
		else {
			Real radian = GrArcCos(
				a.DotProduct(b) *
				a.LengthReciprocal() *
				b.LengthReciprocal()
			);
			if (a.CrossProduct(b) > 0.0)
				return radian;
			else
				return -radian;
		}
	}

	Real VectorAngle(const Vector2& start, const Vector2& end)
	{
		if (start.IsZero() || end.IsZero())
			return 0.0;
		Vector2 a(start);
		Vector2 b(end);
		if (a.LengthSquare() < 1.0) {
			a *= 10000.0;
		}
		if (b.LengthSquare() < 1.0) {
			b *= 10000.0;
		}
		if (IsVectorParallel(a, b)) {
			if (a.DotProduct(b) > 0.0) {
				return 0.0;
			}
			else {
				return 180.0;
			}
		}
		else {
			Real angle = GrArcCosA(
				a.DotProduct(b) *
				a.LengthReciprocal() *
				b.LengthReciprocal()
			);
			if (a.CrossProduct(b) > 0.0)
				return angle;
			else
				return -angle;
		}
	}

} //namespace grid