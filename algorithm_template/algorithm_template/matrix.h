//@header : matrix.h
//@namespace : alg
//@author : zhong
//@date : 2018/9/26

//TODO(zhong) : 此矩阵尚不完善 还有很多矩阵相关的功能没有实现
//				仅仅为了实现矩阵链而实现此矩阵
//				以后看心情补全

#ifndef _ALGORITHM_TEMPLATE_MATRIX_H_
#define _ALGORITHM_TEMPLATE_MATRIX_H_

#include<cassert>
#include<initializer_list>
#include<iostream>

namespace alg
{
	template<typename _Ty>
	class matrix
	{
	public:
		using value_type = _Ty;
		using size_type = size_t;
		using difference_type = ptrdiff_t;
		using pointer = _Ty*;
		using const_pointer = const _Ty*;
		using reference = _Ty& ;
		using const_reference = const _Ty&;

		using iterator = pointer;
		using const_iterator = const_pointer;

		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	public:
		matrix() :	
			m_(nullptr), row_(0), col_(0),size_(0)
		{	}

		matrix(size_type row, size_type col) :
			m_(nullptr), row_(row), col_(col), size_(row * col)
		{
			assert(size_ > 0u);
			assert(m_ = Alloc_(size_, _Ty()));
		}

		matrix(size_type row, size_type col, const _Ty& val) :
			m_(nullptr), row_(row), col_(col), size_(row * col)
		{
			assert(size_ > 0u);
			assert(m_ = Alloc_(size_, val));
		}

		matrix(size_type row, size_type col, std::initializer_list<_Ty> il) :
			m_(nullptr), row_(row), col_(col), size_(row * col)
		{
			assert(size_ > 0u);
			assert(il.size() == size_);
			assert(m_ = Alloc_(size_, _Ty()));
			pointer ptr = m_;
			for (auto& e : il)
			{
				*ptr++ = e;
			}
		}

		matrix(const matrix<_Ty>& right) :
			m_(nullptr), row_(0), col_(0), size_(0)
		{
			if (right.size_ > 0) {
				assert(m_ = Alloc_(right.size_));
				row_ = right.row_;
				col_ = right.col_;
				size_ = right.size_;
				for (size_type i = 0; i < size_; ++i)
				{
					m_[i] = right.m_[i];
				}
			}		
		}

		matrix(matrix<_Ty>&& right) :
			m_(right.m_), row_(right.row_), 
			col_(right.col_), size_(right.size_)
		{
			right.m_ = nullptr;
			right.row_ = 0;
			right.col_ = 0;
			right.size_ = 0;
		}

		matrix<_Ty>& operator=(const matrix<_Ty>& right)
		{
			if (this == &right)
				return *this;
			m_ = reAlloc_(right.size_);
			row_ = right.row_;
			col_ = right.col_;
			size_ = right.size_;
			for (size_type i = 0; i < size; ++i)
			{
				m_[i] = right.m_[i];
			}
			return *this;
		}

		matrix<_Ty>& operator=(matrix<_Ty>&& right)
		{
			if (this == &right)
				return *this;
			Destruct_();
			m_ = right.m_;		
			row_ = right.row_;
			col_ = right.col_;
			size_ = right.size_;
			right.m_ = nullptr;
			right.row_ = 0;
			right.col_ = 0;
			right.size_ = 0;
			return *this;
		}

		matrix<_Ty>& operator=(std::initializer_list<_Ty> il)
		{
			assert(il.size() == size_);
			pointer ptr = m_;
			for (auto& e : il)
			{
				*ptr++ = e;
			}
			return *this;
		}

		void assgin(const _Ty& value)
		{	//assgin value to all elements
			for (size_t i = 0; i < size_; ++i)
			{
				m_[i] = value;
			}
		}

		void fill(const _Ty& value)
		{	// fill matrix by the value
			assgin(value);
		}

		inline iterator begin() noexcept
		{
			return m_;
		}

		inline const_iterator begin() const noexcept
		{
			return m_;
		}

		inline iterator end() noexcept
		{
			return m_ + size_;
		}

		inline const_iterator end() const noexcept
		{
			return m_ + size_;
		}

		inline reverse_iterator rbegin() noexcept
		{
			return reverse_iterator(end());
		}

		inline const_reverse_iterator rbegin() const noexcept
		{
			return const_reverse_iterator(end());
		}

		inline reverse_iterator rend() noexcept
		{
			return reverse_iterator(begin());
		}

		inline const_reverse_iterator rend() const noexcept
		{
			return const_reverse_iterator(begin());
		}

		inline const_iterator cbegin() const noexcept
		{
			return begin();
		}

		inline const_iterator cend() const noexcept
		{
			return end();
		}

		inline const_reverse_iterator crbegin() const noexcept
		{
			return rbegin();
		}

		inline const_reverse_iterator crend() const noexcept
		{
			return rend();
		}

		inline pointer unchecked_begin() noexcept
		{
			return m_;
		}

		inline const_pointer unchecked_begin() const noexcept
		{
			return m_;
		}

		inline pointer unchecked_end() noexcept
		{
			return m_ + size_;
		}

		inline const_pointer unchecked_end() const noexcept
		{
			return m_ + size_;
		}

		inline constexpr size_type size() const noexcept
		{
			return size_;
		}

		inline constexpr size_type max_size() const noexcept
		{
			return size_;
		}

		inline constexpr bool empty() const noexcept
		{	//if m_ is nullptr , then the matrix is empty
			return !m_;
		}

		inline pointer operator[](size_type i) const noexcept
		{
			return m_ + col_ * i;
		}

		bool operator==(const matrix<_Ty>& right) const
		{
			if (row_ == right.row_ && col_ == right.col_){
				bool flag = true;
				for (size_t i = 0; i < size_; ++i)
				{
					if (m_[i] != right.m_[i]){
						flag = false;
						break;
					}
				}
				return flag;
			}
			else {
				return false;
			}
		}

		bool operator!=(const matrix<_Ty>& right) const
		{
			return !(*this == right);
		}

		inline constexpr const_reference at(size_type i, size_type j) const
		{
			assert(i * col_ + j < size_);
			return m_[i * col_ + j];
		}

		inline constexpr const_reference at(size_type i) const
		{
			static_assert(i < size_);
			return m_[i];
		}

		inline reference at(size_type i, size_type j)
		{
			static_assert(i * col_ + j < size_);
			return m_[i * col_ + j];
		}

		inline reference at(size_type i)
		{
			static_assert(i < size_);
			return m_[i];
		}

		inline constexpr size_type row() const noexcept
		{
			return row_;
		}

		inline constexpr size_type column() const noexcept
		{
			return col_;
		}

		inline friend void swap(matrix& left, matrix& right)
		{
			std::swap(left.m_, right.m_);
			std::swap(left.row_, right.row_);
			std::swap(left.col_, right.col_);
		}

		inline reference front() noexcept
		{
			static_assert(m_);
			return m_[0];
		}

		inline const_reference front() const noexcept
		{
			static_assert(m_);
			return m_[0];
		}

		inline reference back() noexcept
		{
			static_assert(m_ & size_ > 0);
			return m_[size_ - 1];
		}

		inline const_reference back() const noexcept
		{
			static_assert(m_ & size_ > 0);
			return m_[size_ - 1];
		}

		inline pointer data() noexcept
		{
			return m_;
		}

		inline const_pointer data() const noexcept
		{
			return m_;
		}

		matrix<_Ty> operator+(const matrix<_Ty>& right) const
		{
			able_to_add(*this, right);
			matrix<_Ty> result(row_, col_);
			for (size_type i = 0; i < size_; ++i)
			{
				result.m_[i] = m_[i] + right.m_[i];
			}
			return result;
		}

		matrix<_Ty> operator-(const matrix<_Ty>& right) const
		{
			able_to_add(*this, right);
			matrix<_Ty> result(row_, col_);
			for (size_type i = 0; i < size_; ++i)
			{
				result.m_[i] = m_[i] - right.m_[i];
			}
			return result;
		}

		matrix<_Ty> operator*(double k) const
		{
			matrix<_Ty> result(*this);
			for (size_type i = 0; i < size_; ++i)
			{
				m_[i] *= k;
			}
			return result;
		}

		matrix<_Ty> operator*(const matrix<_Ty>& right) const
		{
#if _DEBUG
			able_to_multi(*this, right);
#endif // _DEBUG			
			matrix<_Ty> result(row_, right.col_, 0);
			for (size_type i = 0; i < row_; ++i)
			{
				for (size_type j = 0; j < right.col_; ++j)
				{
					for (size_type k = 0; k < col_; ++k)
					{
						result[i][j] += m_[i * col_ + k] * right[k][j];
					}
				}
			}
			return result;
		}

		matrix<_Ty>& operator+=(const matrix<_Ty>& right)
		{
			able_to_add(*this, right);
			for (size_type i = 0; i < size_; ++i)
			{
				m_[i] += right.m_[i];
			}
			return *this;
		}

		matrix<_Ty>& operator-=(const matrix<_Ty>& right)
		{
			able_to_add(*this, right);
			for (size_type i = 0; i < size_; ++i)
			{
				m_[i] -= right.m_[i];
			}
			return *this;
		}

		matrix<_Ty>& operator*=(double k)
		{
			for (size_type i = 0; i < size_; ++i)
			{
				m_[i] *= k;
			}
			return *this;
		}

		friend std::ostream& operator<<(std::ostream& os, const matrix& right) noexcept
		{
			for (size_t i = 0; i < right.row(); ++i)
			{
				os << "[";
				for (size_t j = 0; j < right.column(); ++j)
				{
					os << (right[i][j]) << ((j == right.column() - 1) ? "" : ",");
				}
				os << "]" << std::endl;
			}
			return os;		
		}

	private:

		_Ty * Alloc_(size_type size, const _Ty& val = _Ty())
		{
			return new _Ty[size]{ val };
		}

		_Ty* reAlloc_(size_type size, const _Ty& val = _Ty())
		{
			if (size != size_){
				Destruct_();
				return Alloc_(size, val);
			}
			else {
				return m_;
			}
		}

		void Destruct_()
		{	
			if (m_) {
				size_ == 1u ? delete m_ : delete[] m_;
				row_ = 0;
				col_ = 0;
				size_ = 0;
			}			
		}

		[[noreturn]] void able_to_add(
			const matrix<_Ty>& left, const matrix<_Ty>& right
		) const noexcept
		{
			static_assert(left.row_ == right.row_ && 
				left.col_ == right.col_,
				"tow matrixs can not calculating together");
		}

		[[noreturn]] void able_to_multi(
			const matrix<_Ty>& left, const matrix<_Ty>& right
		) const noexcept
		{
			assert(left.column() == right.row()
				//,"tow matrixs can not calculating together"
			);
		}

	private:

		_Ty*		m_;
		size_type	row_;
		size_type	col_;
		size_type	size_;
	}; //class matrix<_Ty>



} //namespace alg

#endif //_ALGORITHM_TEMPLATE_MATRIX_H_