//@header : immutable_matrix.h
//@author : zhong
//@date : 2018/9/29

#ifndef _ALGORITHM_TEMPLATE_IMMUTABLE_MATRIX_H_
#define _ALGORITHM_TEMPLATE_IMMUTABLE_MATRIX_H_

#include<iostream>
#include<initializer_list>

namespace alg
{
	using real = double;

	template<typename _Ty,
		size_t _row, size_t _col>
	class Matrix
	{
	public:

		using value_type = _Ty;
		using size_type = size_t;

		Matrix() : m{ 0 } {  }

		Matrix(std::initializer_list<_Ty> il) : m{ 0 }
		{
			list_initialization(il);
		}

		Matrix<_Ty, _row, _col>& operator=(std::initializer_list<_Ty> il)
		{
			list_initialization(il);
			return *this;
		}

		constexpr size_t row() const noexcept
		{
			return _row;
		}

		constexpr size_t column() const noexcept
		{
			return _col;
		}

		constexpr size_t size() const noexcept
		{
			return _row * _col;
		}

		_Ty* operator[](int index)
		{
			return m[index];
		}

		_Ty at(int row, int col) const
		{
			checkd_index(row, col);
			return m[row][col];
		}

		_Ty& at(int row, int col)
		{
			checkd_index(row, col);
			return m[row][col];
		}

		_Ty* begin()
		{
			return &m[0][0];
		}

		_Ty* end()
		{
			return begin() + _row * _col;
		}

		const _Ty* begin() const
		{
			return &m[0][0];
		}

		const _Ty* end() const
		{
			return begin() + _row * _col;
		}

		Matrix<_Ty, _row, _col> operator+(const Matrix<_Ty, _row, _col>& right) const
		{
			Matrix<_Ty, _row, _col> result;
			for (size_t i = 0; i < _row; ++i)
			{
				for (size_t j = 0; j < _col; ++j)
				{
					result.m[i][j] = m[i][j] + right.m[i][j];
				}
			}
			return result;
		}

		Matrix<_Ty, _row, _col> operator-(const Matrix<_Ty, _row, _col>& right) const
		{
			Matrix<_Ty, _row, _col> result;
			for (size_t i = 0; i < _row; ++i)
			{
				for (size_t j = 0; j < _col; ++j)
				{
					result.m[i][j] = m[i][j] - right.m[i][j];
				}
			}
			return result;
		}

		Matrix<_Ty, _row, _col> operator*(real k) const
		{
			Matrix<_Ty, _row, _col> result;
			for (size_t i = 0; i < _row; ++i)
			{
				for (size_t j = 0; j < _col; ++j)
				{
					result.m[i][j] = m[i][j] * k;
				}
			}
			return result;
		}

		template<size_t r_col>
		Matrix<_Ty, _row, r_col> operator*(const Matrix<_Ty, _col, r_col>& right) const
		{
			Matrix<_Ty, _row, r_col> result;
			for (size_t i = 0; i < _row; ++i)
			{
				for (size_t j = 0; j < r_col; ++j)
				{
					for (size_t k = 0; k < _col; ++k)
					{
						result.m[i][j] += m[i][k] * right.m[k][j];
					}
				}
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix<_Ty, _row, _col>& matrix) noexcept
		{
			for (size_t i = 0; i < matrix.row(); ++i)
			{
				os << "[";
				for (size_t j = 0; j < matrix.column(); ++j)
				{
					os << matrix.m[i][j] << ((j == matrix.column() - 1) ? "" : ",");
				}
				os << "]" << std::endl;
			}
			return os;
		}

	private:

		void list_initialization(std::initializer_list<_Ty> il)
		{
			assert(_row * _col == il.size());
			_Ty* arr = &m[0][0];
			for (size_t i = 0; i < il.size(); ++i)
			{
				arr[i] = *(il.begin() + i);
			}
		}

		void checkd_index(size_t row, size_t col)
		{
			assert(row < _row);
			assert(col < _col);
		}

		_Ty m[_row][_col];

	}; //class Matrix<_Ty,_row,_col>

} //namespace alg

#endif //_ALGORITHM_TEMPLATE_IMMUTABLE_MATRIX_H_
