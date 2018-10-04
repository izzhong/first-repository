//@header	:	gr_matrix.h
//@abstract	:	1.Define @class Matrix2 and @class Matrix3.
//				2.Implement some math @function relate to Matrix2 or Matrix3 transformation
//@author	:	zhong
//@date		:	2018/8
//@version	:	bate1.0

//@add		:	void CovarianceMatrix2(int size, Vector2* vector_array,Matrix2* covariance_matrix);
//				void CovarianceMatrix3(int size, Vector3* vector_array,Matrix3* covariance_matrix);

//TODO(zhong) : 标准正交基 施密特正交化
//				用于求解向量空间的标准正交基 
//				利用标准正交基可以求出向量在此正交基之下的坐标
//				如果以后用到了可以实现出来 
//				线性代数 p116

#ifndef _GRID_MATRIX_H_
#define _GRID_MATRIX_H_

#include"gr_vector.h"

namespace grid
{
	//@class-->

	//@abstract :	use it like 2d-matrix
	//					m00	 m01
	//					m10	 m11
	//@tips     :	1.use .Transform to transform 2d-vector.
	//				2.use matrix<<=m0,m1,m2,m3; to set a 2d-matrix.
	class Matrix2
	{
	public:
		// m0 m1
		// m2 m3
		Real m[4];
	public:
		//@abstract : 
		//default 2d-matrix is : E
		Matrix2() : 
			m{	1.0,0.0,
				0.0,1.0  } 
		{ }

		explicit Matrix2(
			Real m00, Real m01,
			Real m10, Real m11
		) :
			m{ m00,m01,m10,m11 }
		{ }

		//@abstract : set a 2d-matrix by appoint every component
		Matrix2& Set(Real m00, Real m01, Real m10, Real m11)
		{
			m[0] = m00;
			m[1] = m01;
			m[2] = m10;
			m[3] = m11;
			return *this;
		}

		//@abstract : set a 2d-matrix by another matrix
		Matrix2& Set(const Matrix2& matrix)
		{
			return *this = matrix;
		}

		//@abstract : clear a 2d matrix to zero matrix
		Matrix2& Clear()
		{
			m[0] = m[1] = m[2] = m[3] = 0.0;
			return *this;
		}

		//@abstract : calculate the determinant of the 2d-matrix
		Real Determinant() const
		{
			return m[0] * m[3] - m[1] * m[2];
		}

		//@abstract : check the 2d-matrix is invertible or not
		bool IsInvertible() const
		{	//remove the precision error
			return !IsZeroReal(Determinant());
		}

		//@abstract : invert the 2d-matrix
		//@tips		: uninvertable 2d-matrix will assert
		Matrix2& Invert()
		{
			Real det = Determinant();
			assert(!IsZeroReal(det));
			det = 1.0 / det;
			Matrix2 temp(
				m[3], -m[1],
				-m[2], m[0]
			);
			return *this = temp * det;
		}

		//@abstract : transpose the 2d-matrix
		Matrix2& Transpose()
		{
			Real temp = m[1];
			m[1] = m[2];
			m[2] = temp;
			return *this;
		}

		//@abstract : transform a 2d-vector
		//@para		: the pointer of a 2d-vector which need to transform
		//@tips		: nullptr will assert
		void Transform(Vector2* vector) const
		{
			assert(vector);
			*vector = *this * (*vector);
		}

		//@abstract : transform a 2d-vector array
		//@para		: the pointer to the first 2d-vector.
		//			  the count of the 2c-vectors	
		//@tips		: nullptr will assert
		void Transform(Vector2* vector_array, int size) const
		{
			for (int i = 0; i < size; ++i)
			{
				Transform(&vector_array[i]);
			}
		}

		//@abstract : only calculate real engin value.
		//			  if the matrix have not real engin value
		//			  function will return false
		bool EnginValue(Real* enginvalue_1, Real* enginvalue_2) const
		{
			//assert(enginvalue_1);
			//assert(enginvalue_2);
			Real delta = _real((m[0] + m[3])*(m[0] + m[3]) - 4.0 * Determinant());
			if (delta < 0.0) {
				return false;
			}
			else if (delta == 0.0) {
				if (enginvalue_1 && enginvalue_2) {
					*enginvalue_1 = *enginvalue_2 = (m[0] + m[3])*0.5;
				}
				return true;
			}
			else {
				if (enginvalue_1 && enginvalue_2) {
					*enginvalue_1 = (m[0] + m[3] - GrSqrt(delta)) * 0.5;
					*enginvalue_2 = (m[0] + m[3] + GrSqrt(delta)) * 0.5;
				}
				return true;
			}
		}

		bool EnginVector(Vector2* enginvector_1, Vector2* enginvector_2) const
		{
			assert(enginvector_1);
			assert(enginvector_2);
			Real enginvalue[2] = { 0.0,0.0 };
			Vector2* enginvector[2] = { enginvector_1,enginvector_2 };
			Matrix2 matrix;
			if (EnginValue(&enginvalue[0], &enginvalue[1])) {
				for (int i = 0; i < 2; ++i)
				{
					matrix.Set(
						m[0] - enginvalue[i], m[1],
						m[2], m[3] - enginvalue[i]
					);
					if (!IsZeroReal(matrix.Determinant())) {
						return false;
					}
					enginvector[i]->x = -matrix.m[1];
					enginvector[i]->y = matrix.m[0];
					if (enginvector_1->IsZero()) {
						return false;
					}
					else {
						enginvector[i]->Normalize();
					}
				}
			}
			else {
				return false;
			}
			enginvector_2->Invert();
			return true;
		}

		//@operator overloading-->

		//@abstract	: @operator,
		//			  two funtions use together.
		//			  matrix<<=m0,m1,m2,m3;
		//@warning  : must use it like this.overflow will assert.
		Matrix2& operator<<=(Real r)
		{
			index_ = 0;
			m[index_] = r;
			return *this;
		}

		//@operator<<=
		Matrix2& operator,(Real r)
		{			
			assert(++index_ < 4);
			m[index_] = r;
			return *this;
		}

		//@abstract : 2d-matrix addition
		Matrix2 operator+(const Matrix2& matrix)const
		{
			return Matrix2(
				m[0] + matrix.m[0],
				m[1] + matrix.m[1],
				m[2] + matrix.m[2],
				m[3] + matrix.m[3]
			);
		}

		//@abstract : 2d-matrix subtraction
		Matrix2 operator-(const Matrix2& matrix)const
		{
			return Matrix2(
				m[0] - matrix.m[0],
				m[1] - matrix.m[1],
				m[2] - matrix.m[2],
				m[3] - matrix.m[3]
			);
		}

		//@abstract : negative 2d-matrix
		Matrix2 operator-() const
		{
			return Matrix2(
				-m[0], -m[1],
				-m[2], -m[3]
			);
		}

		//@abstract : 2d-matrix multiplication
		Matrix2 operator*(const Matrix2& matrix)const
		{
			return Matrix2(
				m[0] * matrix.m[0] + m[1] * matrix.m[2],
				m[0] * matrix.m[1] + m[1] * matrix.m[3],
				m[2] * matrix.m[0] + m[3] * matrix.m[2],
				m[2] * matrix.m[1] + m[3] * matrix.m[3]
			);
		}

		//@abstract : 2d-matrix scalar-multiplication
		Matrix2 operator*(Real r)const
		{
			return Matrix2(
				r*m[0], r * m[1], r * m[2], r * m[3]
			);
		}

		//@abstract : 2d-matrix multiplicate 2d-vector
		Vector2 operator*(const Vector2& vector) const
		{
			return Vector2(
				m[0] * vector.x + m[1] * vector.y,
				m[2] * vector.x + m[3] * vector.y
			);
		}

		Matrix2& operator+=(const Matrix2& matrix)
		{
			m[0] += matrix.m[0];
			m[1] += matrix.m[1];
			m[2] += matrix.m[2];
			m[3] += matrix.m[3];
			return *this;
		}

		Matrix2& operator-=(const Matrix2& matrix)
		{
			m[0] -= matrix.m[0];
			m[1] -= matrix.m[1];
			m[2] -= matrix.m[2];
			m[3] -= matrix.m[3];
			return *this;
		}

		Matrix2& operator*=(const Matrix2& matrix)
		{
			Matrix2 temp((*this) * matrix);
			return *this = temp;
		}

		Matrix2& operator*=(Real r)
		{
			m[0] *= r;
			m[1] *= r;
			m[2] *= r;
			m[3] *= r;
			return *this;
		}

		//@abstract : compare two 2d-matrix is equal or not
		bool operator==(const Matrix2& matrix) const
		{
			return
				IsZeroReal(m[0] - matrix.m[0]) &&
				IsZeroReal(m[1] - matrix.m[1]) &&
				IsZeroReal(m[2] - matrix.m[2]) &&
				IsZeroReal(m[3] - matrix.m[3]);
		}

		bool operator!=(const Matrix2& matrix) const
		{
			return !(*this == matrix);
		}

		bool IsZero() const
		{
			return
				IsZeroReal(m[0]) &&
				IsZeroReal(m[1]) &&
				IsZeroReal(m[2]) &&
				IsZeroReal(m[3]);
		}

		bool IsSymmetrical() const
		{
			Matrix2 transpose(*this);
			transpose.Transpose();
			return transpose == *this;
		}

		bool IsDiagonal() const
		{
			return IsZeroReal(m[1]) && IsZeroReal(m[2]);
		}

		bool IsUnit() const
		{
			return
				IsDiagonal() &&
				IsZeroReal(m[0] - 1.0) &&
				IsZeroReal(m[3] - 1.0);
		}

		bool IsOrthogonal() const
		{
			Matrix2 ori(*this);
			ori.Transpose();
			ori *= *this;
			return ori.IsUnit();
		}
	private:
		//@abstract : support the operator<<= and operator,
		int index_ = 0;
	}; //class Matrix2

	//@abstract :	use it like 3d-matrix
	//					m00	 m01  m02
	//					m10	 m11  m12
	//					m20  m21  m22
	//@tips     :	1.use .Transform to transform 3d-vector.
	//				2.use matrix<<=m0,m1,m2,m3,...,m8; to set a 3d-matrix.
	class Matrix3
	{
	public:
		// m0 m1 m2
		// m3 m4 m5
		// m6 m7 m8
		Real m[9];
	public:
		//@abstract : default 3d-matrix is E
		Matrix3() :
			m{  1.0,0.0,0.0,
			    0.0,1.0,0.0,
		        0.0,0.0,1.0  }
		{	}

		explicit Matrix3(
			Real m00, Real m01, Real m02,
			Real m10, Real m11, Real m12,
			Real m20, Real m21, Real m22
		) : 
			m{	m00,m01,m02,
				m10,m11,m12,
				m20,m21,m22  }
		{	}

		//@abstract : set a 3d-matrix by appoint every component
		Matrix3& Set(
			Real m00, Real m01, Real m02,
			Real m10, Real m11, Real m12,
			Real m20, Real m21, Real m22
		){
			m[0] = m00; m[1] = m01; m[2] = m02;
			m[3] = m10; m[4] = m11; m[5] = m12;
			m[6] = m20; m[7] = m21; m[8] = m22;
			return *this;
		}

		//@abstract : set a 3d-matrix by another 3d-matrix
		Matrix3& Set(const Matrix3& matrix)
		{
			*this = matrix;
		}

		//@abstract : clear a 3d matrix to zero matrix
		Matrix3& Clear()
		{
			m[0] = m[1] = m[2] = 
			m[3] = m[4] = m[5] = 
			m[6] = m[7] = m[8] = 0;
			return *this;
		}

		//@abstract : calculate the determinant of a 3d-matrix
		Real Determinant() const
		{
			return
				(m[3] * m[7] - m[4] * m[6]) * m[2] +
				(m[5] * m[6] - m[3] * m[8]) * m[1] +
				(m[4] * m[8] - m[5] * m[7]) * m[0];
		}

		//@abstract : check a 3d-matrix is invertible or not
		bool IsInvertible() const
		{
			return !IsZeroReal(Determinant());
		}

		//@abstract : invert a 3d-matrix
		//@tips		: uninvertible matrix will assert
		Matrix3& Invert()
		{
			Real det = Determinant();
			assert(!IsZeroReal(det));
			det = 1.0 / det;
			Matrix3 temp(
				m[4] * m[8] - m[5] * m[7],
				m[2] * m[7] - m[1] * m[8],
				m[1] * m[5] - m[2] * m[4],
				m[5] * m[6] - m[3] * m[8],
				m[0] * m[8] - m[2] * m[6],
				m[2] * m[3] - m[0] * m[5],
				m[3] * m[7] - m[4] * m[6],
				m[1] * m[6] - m[0] * m[7],
				m[0] * m[4] - m[1] * m[3]
			);
			return *this = temp * det;
		}

		bool EnginValue(Real* enginvalue_1, Real* enginvalue_2) const
		{
			Real delta = (m[0] - m[3])*(m[0] - m[3]) - 4.0*m[1] * m[2];
			if (delta < 0.0) {
				return false;
			}
			else if (IsZeroReal(delta)) {
				*enginvalue_1 = *enginvalue_2 = (m[0] + m[3])*0.5;
				return true;
			}
			else {
				*enginvalue_1 = (m[0] + m[3] + GrSqrt(delta)) * 0.5;
				*enginvalue_2 = (m[0] + m[3] - GrSqrt(delta)) * 0.5;
				return true;
			}
		}

		bool EnginVector(Vector2* enginvector_1, Vector2* enginvector_2) const
		{
			Real enginvalue[2] = { 0.0,0.0 };
			Vector2* enginvector[2] = { enginvector_1,enginvector_2 };
			Matrix2 matrix;
			if (EnginValue(&enginvalue[0], &enginvalue[1])) {
				for (int i = 0; i < 2; ++i) {
					matrix.Set(
						m[0] - enginvalue[i], m[1],
						m[2], m[3] - enginvalue[i]
					);
					if (matrix.IsInvertible()) {
						*enginvector[i] = matrix.Invert() * (Vector2(0.0, 0.0));
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

		//@abstract : transpose a 3d-matrix
		Matrix3& Transpose()
		{
			Real temp = m[1]; m[1] = m[3]; m[3] = temp;
			temp = m[2]; m[2] = m[6]; m[6] = temp;
			temp = m[5]; m[5] = m[7]; m[7] = temp;
			return *this;
		}

		//@abstract : transform a 3d-vector
		//@para		: the pointer of a 3d-vector which need to transform
		//@tips		: nullptr will assert
		void Transform(Vector3* vector) const
		{
			assert(vector);
			*vector = *this * (*vector);
		}

		//@abstract : transform a 3d-vector array
		//@para		: 1.the count of the 3d-vectors 
		//			  2.the pointer to the first 3d-vector.
		//@tips		: nullptr will assert
		void Transform(int size, Vector3* vector_array) const
		{
			for (int i = 0; i < size; ++i)
			{
				Transform(&vector_array[i]);
			}
		}

		//@abstract : transform a 2d-vector
		//			  althought it looks strange
		//			  but it can simplify the calculate process
		//@para		: the pointer to 2d-vector
		//@tips		: it looks the 2d-vector as 3d-vector : (x,y,1.0)
		//			  to run the 3d-matrix multiplication
		void Transform(Vector2* vector) const
		{
			assert(vector);
			*vector = *this * (*vector);
		}

		void Transform(int size, Vector2* vector_array) const
		{
			for (int i = 0; i < size; ++i)
			{
				Transform(&vector_array[i]);
			}
		}

		//@operator overloading-->

		//@Matrix2.operator<<=
		Matrix3& operator<<=(Real r)
		{
			index_ = 0;
			m[index_] = r;
			return *this;
		}

		//@Matrix2.operator,
		Matrix3& operator,(Real r)
		{
			assert(++index_ < 9);
			m[index_] = r;
			return *this;
		}

		Matrix3 operator+(const Matrix3& matrix) const
		{
			return Matrix3(
				m[0] + matrix.m[0], m[1] + matrix.m[1], m[2] + matrix.m[2],
				m[3] + matrix.m[3], m[4] + matrix.m[4], m[5] + matrix.m[5],
				m[6] + matrix.m[6], m[7] + matrix.m[7], m[8] + matrix.m[8]
			);
		}

		Matrix3 operator-(const Matrix3& matrix) const
		{
			return Matrix3(
				m[0] - matrix.m[0], m[1] - matrix.m[1], m[2] - matrix.m[2],
				m[3] - matrix.m[3], m[4] - matrix.m[4], m[5] - matrix.m[5],
				m[6] - matrix.m[6], m[7] - matrix.m[7], m[8] - matrix.m[8]
			);
		}

		Matrix3 operator-() const
		{
			return Matrix3(
				-m[0], -m[1], -m[2],
				-m[3], -m[4], -m[5],
				-m[6], -m[7], -m[8]
			);
		}

		Matrix3 operator*(const Matrix3& matrix) const
		{
			return Matrix3(
				m[0] * matrix.m[0] + m[1] * matrix.m[3] + m[2] * matrix.m[6],
				m[0] * matrix.m[1] + m[1] * matrix.m[4] + m[2] * matrix.m[7],
				m[0] * matrix.m[2] + m[1] * matrix.m[5] + m[2] * matrix.m[8],

				m[3] * matrix.m[0] + m[4] * matrix.m[3] + m[5] * matrix.m[6],
				m[3] * matrix.m[1] + m[4] * matrix.m[4] + m[5] * matrix.m[7],
				m[3] * matrix.m[2] + m[4] * matrix.m[5] + m[5] * matrix.m[8],

				m[6] * matrix.m[0] + m[7] * matrix.m[3] + m[8] * matrix.m[6],
				m[6] * matrix.m[1] + m[7] * matrix.m[4] + m[8] * matrix.m[7],
				m[6] * matrix.m[2] + m[7] * matrix.m[5] + m[8] * matrix.m[8]
			);
		}

		Matrix3 operator*(Real r) const
		{
			return Matrix3(
				r * m[0], r * m[1], r * m[2],
				r * m[3], r * m[4], r * m[5],
				r * m[6], r * m[7], r * m[8]
			);
		}

		Vector3 operator*(const Vector3& vector) const
		{
			return Vector3(
				m[0] * vector.x + m[1] * vector.y + m[2] * vector.z,
				m[3] * vector.x + m[4] * vector.y + m[5] * vector.z,
				m[6] * vector.x + m[7] * vector.y + m[8] * vector.z
			);
		}

		Matrix3& operator+=(const Matrix3& matrix)
		{
			m[0] += matrix.m[0]; m[1] += matrix.m[1]; m[2] += matrix.m[2];
			m[3] += matrix.m[3]; m[4] += matrix.m[4]; m[5] += matrix.m[5];
			m[6] += matrix.m[6]; m[7] += matrix.m[7]; m[8] += matrix.m[8];
			return *this;
		}

		Matrix3& operator-=(const Matrix3& matrix)
		{
			m[0] -= matrix.m[0]; m[1] -= matrix.m[1]; m[2] -= matrix.m[2];
			m[3] -= matrix.m[3]; m[4] -= matrix.m[4]; m[5] -= matrix.m[5];
			m[6] -= matrix.m[6]; m[7] -= matrix.m[7]; m[8] -= matrix.m[8];
			return *this;
		}

		Matrix3& operator*=(const Matrix3& matrix)
		{
			Matrix3 temp((*this) * matrix);
			return *this = temp;
		}

		Matrix3& operator*=(Real r)
		{
			m[0] *= r;	m[1] *= r;	m[2] *= r;
			m[3] *= r;	m[4] *= r;	m[5] *= r;
			m[6] *= r;	m[7] *= r;	m[8] *= r;
			return *this;
		}

		bool operator==(const Matrix3& matrix) const
		{
			return
				IsZeroReal(m[0] - matrix.m[0]) &&
				IsZeroReal(m[1] - matrix.m[1]) &&
				IsZeroReal(m[2] - matrix.m[2]) &&
				IsZeroReal(m[3] - matrix.m[3]) &&
				IsZeroReal(m[4] - matrix.m[4]) &&
				IsZeroReal(m[5] - matrix.m[5]) &&
				IsZeroReal(m[6] - matrix.m[6]) &&
				IsZeroReal(m[7] - matrix.m[7]) &&
				IsZeroReal(m[8] - matrix.m[8]);
		}

		bool operator!=(const Matrix3& matrix) const
		{
			return !(*this == matrix);
		}

		bool IsDiagonal() const
		{
			return 
				IsZeroReal(m[1]) && 
				IsZeroReal(m[2]) &&
				IsZeroReal(m[3]) && 
				IsZeroReal(m[5]) &&
				IsZeroReal(m[6]) && 
				IsZeroReal(m[7]);
		}

		bool IsUnit() const
		{
			return
				IsDiagonal() &&
				IsZeroReal(m[0] - 1.0) &&
				IsZeroReal(m[4] - 1.0) &&
				IsZeroReal(m[8] - 1.0);
		}

		bool IsSymmetrical() const
		{
			Matrix3 transpose(*this);
			transpose.Transpose();
			return transpose == *this;
		}

		bool IsZero() const
		{
			return
				IsZeroReal(m[0]) &&
				IsZeroReal(m[1]) &&
				IsZeroReal(m[2]) &&
				IsZeroReal(m[3]) &&
				IsZeroReal(m[4]) &&
				IsZeroReal(m[5]) &&
				IsZeroReal(m[6]) &&
				IsZeroReal(m[7]) &&
				IsZeroReal(m[8]);
		}

		bool IsOrthogonal() const
		{
			Matrix3 ori(*this);
			ori.Transpose();
			ori *= *this;
			return ori.IsUnit();
		}

	private:
		//@Matrix2.index_
		int index_ = 0;
	private:
		//@abstract : this funtion is only use for simplify matrix transformation
		//			  it see the 2d-vector as 3d-vector(x,y,1.0).
		Vector2 operator*(const Vector2& vector) const
		{
			return Vector2(
				m[0] * vector.x + m[1] * vector.y + m[2],
				m[3] * vector.x + m[4] * vector.y + m[5]
			);
		}
	}; //class Matrix3

	//@function-->

	//@tips		: when you use the functions below
	//@warning	: the output matrix only can use to transform 2d-vector.
	//			  althought the 3d-vector para is legal
	//			  but the result can not be correct.
	//			  use it like
	//				matrix.Transform(size,&vector2);

	//@abstract : by the input para calculate the translate matrix
	//		      and assgin to the given pointer.
	//@para		: input : 2d-vector which stands for the translation.
	//			: ouput : the 3d-matrix pointer which will get the translate-matrix.
	//@tip:		: input 2d-vector can be any vector.
	//			  output pointer must not be a nullptr,otherwise will assert.
	//			  after get the matrix.you can use it like this
	//				matrix.Transform(size,&2d_vector_array);
	//@warning	: the output matrix only can use to transform 2d-vector.
	//			  althought the 3d-vector para is legal
	//			  but the result can not be correct.
	void TranslateMatrix(
		const Vector2& translate,	//input
		Matrix3* matrix				//output
	);

	//@abstract : by the input para calculate the rotate matrix
	//			  and assgin to the given matrix pointer.
	//@para		: input : 1.the center of the rotatetion
	//					  2.the angle of the rotatetion.
	//					    if the angle > 0,the rotatetion will be unti-clock-wise
	//@tips		: input paras have no limits
	//			  output matrix pointer must not be nullptr otherwise will assert
	//@warning  : @TranslateMatrix
	void RotateMatrix(
		const Vector2& center, Real angle,	//input
		Matrix3* matrix						//output
	);

	//@abstract : by the given center and coefficient calculate the
	//			  scale matrix and assgin to the given matrix pointer
	//@para		: input 1.the center of the scale
	//					2.the scale coefficient of x-axis. x = coe * x0.
	//					3.the scale coefficient of y-axis.
	//@tips		: after the transformation.
	//				x -> coe * x0
	//				y -> coe * y0
	//@warning	: @TranslateMatrix
	void ScaleMatrix(
		const Vector2& center, Real coefficient_x, Real coefficient_y,	//input
		Matrix3* matrix													//output
	);

	//@abstract : by the given center the scale normal and coefficient
	//			  calculate the scale matrix
	//@para		: input 1.the center of the scale
	//					2.the scale will along the axis , which normal vector is [normal]
	//					3.the scale coefficient
	//			  ouput : the 3d-matrix pointer which will get the matrix
	//@tips:	  input para [direction] must be a normalization vector.
	//			  otherwise the result can not be correct.
	//@warning	: @TranslateMatrix
	void ScaleMatrix(
		const Vector2& center, const Vector2& normal, Real coefficient,	//input
		Matrix3* matrix														//output
	);

	//@abstract : by the given center and project direction axis
	//			  calculate the project matrix
	//@para		: input 1. projection center
	//					2. projection direction axis vector
	//			  output :
	//@tips		: input para[direction] must ve a normalization vector.
	//			  after the transformation.the order of vectors will be opposite.
	//@warning	: @TranslateMatrix
	void ProjectMatrix(
		const Vector2& center, const Vector2& direction,	//input
		Matrix3* matrix										//output
	);

	//@abstract : by the para calculate the reflection matrix
	//@para		: input 1.the reflecttion center
	//					2.the reflection direction axis
	//@tips		: the para[direction] must be a normalization vector
	//			  after the transformation.the order of vectors will be opposite.
	//@warning	: @TranslateMatrix
	void ReflectMatrix(
		const Vector2& center, const Vector2& direction,	//input
		Matrix3* matrix										//output
	);

	//@abstract : get the matrix which can shear along x-axis
	//@para		: input 1.shear center
	//					2.shear coefficient
	//@tips		: if coefficient = 0 , no shear
	//			  if coefficient > 0 , shear along positive direction of x-axis
	//			  if coefficient < 0 , shear along negative direction of x-axis
	//@warning	: @TranslateMatrix
	void ShearXMatrix( 
		const Vector2& center, Real coefficient,	//input
		Matrix3* matrix								//output
	);

	//@abstract : get the matrix which can shear along y-axis
	//@para		: input 1.shear center
	//					2.shear coefficient
	//@tips		: if coefficient = 0 , no shear
	//			  if coefficient > 0 , shear along positive direction of y-axis
	//			  if coefficient < 0 , shear along negative direction of y-axis
	//@warning	: @TranslateMatrix
	void ShearYMatrix(
		const Vector2& center, Real coefficient,	//input
		Matrix3* matrix								//output
	);

	//@abstract : Matrix Translation Transformation
	//@para		: 1.the translation of the transformation
	//			  2.vector_array_size and the first address
	//			  3.potional para use to get the inverse matrix.
	//			  if want , pass a matrix pointer.
	//@tips		: vector_array shoule be valid.
	//			  otherwise the result can not be correct
	//@warning	: @TranslateMatrix
	void TranslateTransform(
		const Vector2& translate,			//input
		int size, Vector2* vector_array,	//vector array for transformation
		Matrix3* inverse_matrix = nullptr	//optional para : get the inverse matrix
	);

	//@abstract : @TranslateTransform
	//			  @RotateMatrix
	void RotateTransform(
		const Vector2& center, Real angle,	//input
		int size, Vector2* vector_array,	//vector array for transformation
		Matrix3* inverse_matrix = nullptr	//optinal para : get the inverse matrix
	);

	//@abstract : @TranslateTransform
	//			  @ScaleTransform(const Vector2&,Real,Real)
	//@tips		: only the coefficient is not 0.0
	//			  the validly inverse_matrix can be gotton.
	void ScaleTransform(
		const Vector2& center, Real coefficient_x, Real coefficient_y,	//input
		int size, Vector2* vector_array,								//vector array for transformation
		Matrix3* inverse_matrix = nullptr								//optinal para : get the inverse matrix
	);

	//@abstract : @TranslateTransform
	//			  @ScaleMatrix(const Vector2&,const Vector2,Real)
	//@tips		: only the coefficient is not 0.0
	//			  the validly inverse_matrix can be gotton.
	void ScaleTransform(
		const Vector2& center, const Vector2& normal, Real coefficient,	//input
		int size, Vector2* vector_array,								//vector array for transformation
		Matrix3* inverse_matrix	= nullptr								//optinal para : get the inverse matrix						
	);

	//@abstract : @TranslateTransform
	//			  @ProjectionMatrix
	//@tips		: can not get the inverse matrix
	void ProjectTransform(
		const Vector2& center, const Vector2& normal,	//input
		int size, Vector2* vector_array					//vector array for transformation
	);

	//@abstract : @TranslateTransform
	//			  @ReflectMatrix	
	void ReflectTransform(
		const Vector2& center, const Vector2& normal,	//input
		int size, Vector2* vector_array,				//vector array for transformation
		Matrix3* inverse_matrix = nullptr				//optinal para : get the inverse matrix
	);

	//@abstract : @TranslateTransform
	//			  @ShearXMatrix
	void ShearXTransform(
		const Vector2& center, Real coefficient,	//input
		int size, Vector2* vector_array,			//vector array for transformation
		Matrix3* inverse_matrix = nullptr			//optinal para : get the inverse matrix
	);

	//@abstract : @TranslateTransform
	//			  @ShearYMatrix
	void ShearYTransform(
		const Vector2& center, Real coefficient,	//input
		int size, Vector2* vector_array,			//vector array for transformation
		Matrix3* inverse_matrix = nullptr			//optinal para : get the inverse matrix
	);

	//@abstract : 二维协方差矩阵
	void CovarianceMatrix2(
		int size, Vector2* vector_array,
		Matrix2* covariance_matrix
	);

	//@abstract : 三维协方差矩阵
	void CovarianceMatrix3(
		int size, Vector3* vector_array,
		Matrix3* covariance_matrix
	);

} //namespace grid

#endif //_GRID_MATRIX_H_