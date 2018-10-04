//@source   : gr_matrix.cpp
//@abstract : Implement Translate,Rotate,Scale,Project,Reflect,Shear MatrixTransformation.
//@author	: zhong
//@date		: 2018/8
//@version	: bate1.0

//TODO(zhong) : 添加一些特殊矩阵类型 比如 对称 对角 等等矩阵的判断

#include"include/gr_matrix.h"

namespace grid
{
	//@header:gr_matrix.h[TranslateMatrix]
	void TranslateMatrix(
		const Vector2& translate,	//input
		Matrix3* matrix				//output
	){
		assert(matrix);	//nullpte will assert
		matrix->Set(	
			1.0, 0.0, translate.x,
			0.0, 1.0, translate.y,
			0.0, 0.0, 1.0
		);
	}

	//@header:gr_matrix.h[RotateMatrix]
	void RotateMatrix(
		const Vector2& center, Real angle,	//input
		Matrix3* matrix						//ouput
	) {
		assert(matrix);						//nullptr will assert
		angle = AngleToRadian(angle);		//change to the radian
		matrix->m[0] = GrCos(angle);		//para is radian
		matrix->m[1] = -GrSin(angle);
		matrix->m[2] =						//calculate the translation
			center.x - matrix->m[0] * center.x - matrix->m[1] * center.y;
		matrix->m[3] = -matrix->m[1];
		matrix->m[4] = matrix->m[0];
		matrix->m[5] =						//calculate the translation 
			center.y - matrix->m[3] * center.x - matrix->m[4] * center.y;
		matrix->m[6] = matrix->m[7] = 0.0;
		matrix->m[8] = 1.0;
	}

	//@header:gr_matrix.h[ScaleMatrix]
	void ScaleMatrix(
		const Vector2& center, Real coefficient_x, Real coefficient_y,	//input
		Matrix3* matrix													//output
	) {
		assert(matrix);
		matrix->Set(
			coefficient_x, 0.0, center.x - coefficient_x * center.x,
			0.0, coefficient_y, center.y - coefficient_y * center.y,
			0.0, 0.0, 1.0
		);
	}

	void ScaleMatrix(
		const Vector2& center, const Vector2& normal, Real coefficient,
		Matrix3* matrix
	) {
		assert(matrix);
		matrix->m[0] = 1 + (coefficient - 1) * normal.x * normal.x;
		matrix->m[1] = (coefficient - 1) * normal.x * normal.y;
		matrix->m[2] = center.x - matrix->m[0] * center.x - 
			matrix->m[1] * center.y;
		matrix->m[3] = matrix->m[1];
		matrix->m[4] = 1 + (coefficient - 1) * normal.y * normal.y;
		matrix->m[5] = center.y - matrix->m[3] * center.x - 
			matrix->m[4] * center.y;
		matrix->m[6] = matrix->m[7] = 0.0;
		matrix->m[8] = 1.0;
	}

	void ProjectMatrix(
		const Vector2& center, const Vector2& normal,
		Matrix3* matrix
	) {
		ScaleMatrix(center, normal, 0.0, matrix);
	}

	void ReflectMatrix(
		const Vector2& center, const Vector2& normal,
		Matrix3* matrix
	) {
		ScaleMatrix(center, VectorVerticalProduct(normal), -1.0, matrix);
	}

	void ShearYMatrix(
		const Vector2& center, Real coefficient,
		Matrix3* matrix
	) {
		assert(matrix);
		matrix->Set(
			1.0, 0.0, 0.0,
			coefficient, 1.0, -coefficient * center.x,
			0.0, 0.0, 1.0
		);
	}

	void ShearXMatrix(
		const Vector2& center, Real coefficient,
		Matrix3* matrix
	) {
		assert(matrix);
		matrix->Set(
			1.0, coefficient, -coefficient * center.y,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0
		);
	}

	void TranslateTransform(
		const Vector2& translate,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix(
			1.0, 0.0, translate.x,
			0.0, 1.0, translate.y,
			0.0, 0.0, 1.0
		);
		matrix.Transform(size, vector_array);
		if (inverse_matrix)	{
			matrix.m[2] = -translate.x;
			matrix.m[5] = -translate.y;
			*inverse_matrix = matrix;
		}		
	}

	void RotateTransform(
		const Vector2& center, Real angle,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		RotateMatrix(center, angle, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix)
		{
			*inverse_matrix = matrix;
			inverse_matrix->m[1] = -matrix.m[1];
			inverse_matrix->m[2] = center.x - inverse_matrix->m[0] * center.x - 
				inverse_matrix->m[1] * center.y;
			inverse_matrix->m[3] = -matrix.m[3];
			inverse_matrix->m[5] = center.y - inverse_matrix->m[3] * center.x - 
				inverse_matrix->m[4] * center.y;
		}
	}

	void ScaleTransform(
		const Vector2& center, Real coefficient_x, Real coefficient_y,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		ScaleMatrix(center, coefficient_x, coefficient_y, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix && coefficient_x != 0.0 && coefficient_y != 0.0)
		{
			ScaleMatrix(
				center, 1.0 / coefficient_x, 1.0 / coefficient_y, 
				inverse_matrix
			);
		}
	}

	void ScaleTransform(
		const Vector2& center, const Vector2& normal, Real coefficient,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		ScaleMatrix(center, normal, coefficient, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix && coefficient != 0.0)
		{
			ScaleMatrix(center, normal, 1.0 / coefficient, inverse_matrix);
		}
	}

	void ProjectTransform(
		const Vector2& center, const Vector2& normal,
		int size, Vector2* vector_array
	) {
		assert(vector_array);
		Matrix3 matrix;
		ProjectMatrix(center, normal, &matrix);
		matrix.Transform(size, vector_array);
	}

	void ReflectTransform(
		const Vector2& center, const Vector2& normal,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		ReflectMatrix(center, normal, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix)
			*inverse_matrix = matrix;
	}

	void ShearXTransform(
		const Vector2& center, Real coefficient,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		ShearXMatrix(center, coefficient, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix)
		{
			ShearXMatrix(center, -coefficient, inverse_matrix);
		}
	}

	void ShearYTransform(
		const Vector2& center, Real coefficient,
		int size, Vector2* vector_array,
		Matrix3* inverse_matrix
	) {
		assert(vector_array);
		Matrix3 matrix;
		ShearYMatrix(center, coefficient, &matrix);
		matrix.Transform(size, vector_array);
		if (inverse_matrix)
		{
			ShearYMatrix(center, -coefficient, inverse_matrix);
		}
	}

	void CovarianceMatrix2(
		int size, Vector2* vector_array,
		Matrix2* covariance_matrix
	) {
		assert(size > 0);
		assert(vector_array);
		assert(covariance_matrix);
		Vector2 variance = GrVariance(size, vector_array);
		Real covariance = Covariance(size, vector_array);
		covariance_matrix->m[0] = variance.x;
		covariance_matrix->m[1] = covariance;
		covariance_matrix->m[2] = covariance;
		covariance_matrix->m[3] = variance.y;
	}

	void CovarianceMatrix3(
		int size,Vector3* vector_array,
		Matrix3* covariance_matrix
	) {
		assert(size > 0);
		assert(vector_array);
		assert(covariance_matrix);
		Vector3 variance = GrVariance<Vector3>(size, vector_array);
		Vector3 covariance = Covariance(size, vector_array);
		covariance_matrix->m[0] = variance.x;
		covariance_matrix->m[1] = covariance.x;
		covariance_matrix->m[2] = covariance.z;
		covariance_matrix->m[3] = covariance.x;
		covariance_matrix->m[4] = variance.y;
		covariance_matrix->m[5] = covariance.y;
		covariance_matrix->m[6] = covariance.z;
		covariance_matrix->m[7] = covariance.y;
		covariance_matrix->m[8] = variance.z;
	}

} //namespace grid