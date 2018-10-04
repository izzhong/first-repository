//@header : gr_debug_cout.h
//@author : zhong
//@abstract : 用于方便的输出测试元素

//@abstract : 函数分为两种类型 
//				扩展cout输出
//				测试输出

#ifndef _GRID_DEBUG_COUT_H_
#define _GRID_DEBUG_COUT_H_

#include<gr_geometry.h>
#include<iostream>

namespace grid
{
	namespace grid_debug
	{
		std::ostream& operator<<(std::ostream& os, const Vector2& right)
		{
			os << "(" << _real(right.x) << "," << _real(right.y) << ")";
			return os;
		}

		std::ostream& operator<<(std::ostream& os, const Vector3& right)
		{
			os << " (" << _real(right.x) << "," << _real(right.y) << "," << _real(right.z) << ")";
			return os;
		}

		std::ostream& operator<<(std::ostream& os, const PointNode& pointnode)
		{
			os
				<< "[ addr : " << &pointnode << " ]" << std::endl
				<< "[ point : " << pointnode.point << " ]" << std::endl
				<< "[ prev : " << pointnode.prev << " ]" << std::endl
				<< "[ next : " << pointnode.next << " ]" << std::endl;
			return os;
		}

		std::ostream& operator<<(std::ostream& os, const PointArray& pointarray)
		{
			for (int i = 0; i < pointarray.Size(); ++i)
			{
				os << pointarray[i] << " ";
				if (i % 5 == 0) {
					os << std::endl;
				}
			}		
		}

		std::ostream& operator<<(std::ostream& os,const PointSet& pointset)
		{
			int i = 0;
			auto work = pointset.HeadPtr();
			do {
				os << work->point << " ";
				work = work->next;
				if (i++ % 5 == 0) {
					os << std::endl;
				}
			} while (work != pointset.HeadPtr());
			return os;
		}

		std::ostream& operator<<(std::ostream& os, const Matrix2& matrix)
		{
			os
				<< "[" << _real(matrix.m[0]) << "," << _real(matrix.m[1]) << "]" << std::endl
				<< "[" << _real(matrix.m[2]) << "," << _real(matrix.m[3]) << "]" << std::endl;
			return os;
		}

		std::ostream& operator<<(std::ostream& os, const Matrix3& matrix)
		{
			os
				<< "[" << _real(matrix.m[0]) << "," << _real(matrix.m[1]) << "," << _real(matrix.m[2]) << "]" << std::endl
				<< "[" << _real(matrix.m[3]) << "," << _real(matrix.m[4]) << "," << _real(matrix.m[5]) << "]" << std::endl
				<< "[" << _real(matrix.m[6]) << "," << _real(matrix.m[7]) << "," << _real(matrix.m[8]) << "]" << std::endl;
			return os;
		}

		template<typename T>
		std::ostream& CoutAddress(std::ostream& os, const T& t)
		{
			os << "[ adr : " << &t << " ] ";
			return os;
		}

		template<typename T>
		std::ostream& Cout(std::ostream& os, const T& t)
		{
			CoutAddress(os, t) << t;
			return os;
		}

		std::ostream& CoutRealArray(std::ostream& os, int size, Real* real_array)
		{
			if (real_array) {
				for (int i = 0; i < size; ++i)
				{
					CoutAddress(std::cout, real_array[i])
						<< real_array[i] << std::endl;
				}
			}
			else {
				os << "[nullptr]" << std::endl;
			}
			return os;
		}

		template<typename Vector>
		std::ostream& CoutVectorArray(std::ostream& os, int size, Vector* vector_array)
		{
			if (vector_array) {
				for (int i = 0; i < size; ++i)
				{
					os << vector_array[i] << std::endl;
				}
			}
			else {
				os << "[nullptr]" << std::endl;
			}
			return os;
		}

		std::ostream& CoutTest(std::ostream& os, const Vector2& right)
		{
			os
				<< " [" << " val : " << " (" << _real(right.x) << "," << _real(right.y) << ")" << " ] "
				<< " [" << " len : " << right.Length() << " ] "
				<< " [" << " len2 : " << right.LengthSquare() << " ] "
				<< " [" << " zero? : " << right.IsZero() << " ] "
				<< " [" << " nor? : " << right.IsNormalization() << " ] ";
			if (!right.IsZero()) {
				os << " [" << " lenr : " << right.LengthReciprocal() << " ] ";
				Vector2 temp(right);
				temp.Normalize();
				os
					<< " [" << " nor : " << " (" << _real(temp.x) << "," << _real(temp.y) << ")" << " ]"
					<< " [" << " norlen : " << temp.Length() << " ] "
					<< std::endl;
			}	
			return os;
		}

		std::ostream& CoutTest(std::ostream& os, const Vector3& right)
		{
			os
				<< " [" << " val : " << " (" << _real(right.x) << "," << _real(right.y) << "," << _real(right.z) << ")" << " ] "
				<< " [" << " len : " << right.Length() << " ] "
				<< " [" << " len2 : " << right.LengthSquare() << " ] "
				<< " [" << " zero? : " << right.IsZero() << " ] "
				<< " [" << " nor? : " << right.IsNormalization() << " ] ";
			if (!right.IsZero()) {
				os << " [" << " lenr : " << right.LengthReciprocal() << " ] ";
				Vector3 temp(right);
				temp.Normalize();
				os
					<< " [" << " nor : " << " (" << _real(temp.x) << "," << _real(temp.y) << "," << _real(temp.z) << ")" << " ]"
					<< " [" << " norlen : " << temp.Length() << " ] "
					<< std::endl;
			}
			return os;
		}

		template<typename Matrix>
		std::ostream& CoutTest(std::ostream& os , const Matrix& matrix)
		{
			Matrix m(matrix);
			CoutAddress(os, matrix) << endl
				<< matrix
				<< "[ det : " << matrix.Determinant() << " ] " << std::endl
				<< "[ transpose ]" << std::endl
				<< m.Transpose()
				<< "[ invertible? : " << matrix.IsInvertible() << " ] " << std::endl;
			if (matrix.IsInvertible()) {
				m.Transpose().Invert();
				os
					<< "[ invert ]" << std::endl
					<< m;
			}
			return os;
		}

	} //namespace grid_debug

} //namespace grid

#endif //_GRID_DEBUG_COUT_H_
