//@header : matrix_chain_multiplication.h
//@author : zhong
//@date : 2018/9/30

//TODO(zhong) : 
//		目前已经想到的两个优化方向:
//		1.过程矩阵 S M 的结构使得其仅需要半矩阵就可存储 而程序为了方便使用了全矩阵
//		2.矩阵链乘法的计算使用递归实现，资源占用量超大

#ifndef _ALGORITHM_TEMPLATE_MATRIX_CHAIN_MULTIPLICATION_H_
#define _ALGORITHM_TEMPLATE_MATRIX_CHAIN_MULTIPLICATION_H_

#include"matrix.h"

namespace alg
{
	using size_s = size_t;
	using size_m = long long;

	template<typename _Matrixl, typename _Matrixr>
	auto MatrixMultiplicationComplexity(const _Matrixl& left, const _Matrixr& right) noexcept -> typename _Matrixl::size_type
	{
		return left.row() * left.column() * right.column();
	}

	template<typename _Matrix_Chain>
	long long MatrixMultiplicationComplexity(const _Matrix_Chain& matrix_chain)
	{
		long long complexity = 0;
		size_t size = matrix_chain.size();
		if (size == 0 || size == 1) {
			return 0;
		}
		else {
			for (size_t i = 0; i < size - 1; ++i)
			{
				complexity += matrix_chain.begin()->row() * (matrix_chain.begin() + i)->column() * (matrix_chain.begin() + i + 1)->column();
			}
			return complexity;
		}
	}

	template<typename _Matrix_Chain, typename _Matrix>
	std::pair<size_s, size_m> MatrixChainPartition(const _Matrix_Chain& matrix_chain, const _Matrix& M, size_t i, size_t j)
	{
		assert(i <= j && j < matrix_chain.size());
		if (i == j) {
			return { i,0 };
		}
		else if (j - i == 1) {
			return { i,MatrixMultiplicationComplexity(*(matrix_chain.begin() + i),*(matrix_chain.begin() + j))};
		}
		else { 
			size_s s = i;
			size_m m = M[i][i] + M[i + 1][j] + MatrixMultiplicationComplexity(*(matrix_chain.begin() + i), *(matrix_chain.begin() + j));
			size_m t = 0;
			for (size_t k = i + 1; k < j; ++k)
			{
				t = M[i][k] + M[k + 1][j] + MatrixMultiplicationComplexity(*(matrix_chain.begin() + i), *(matrix_chain.begin() + j));
				if (t < m) {
					s = k;
					m = t;
				}
			}
			return { s,m };
		}
	}

	template<typename _Matrix_Chain>
	auto MatrixChainMultiplication(const _Matrix_Chain& matrix_chain, const matrix<size_t>& S, size_t i, size_t j) -> typename _Matrix_Chain::value_type
	{
		assert(i <= j);
		if (i == j) {
			return *matrix_chain.begin();
		}
		if (j - i == 1) {
			return *matrix_chain.begin() * *(matrix_chain.begin() + 1);
		}
		size_s partition = S.at(i, j);
		_Matrix_Chain MC1(matrix_chain.begin(), matrix_chain.begin() + (partition - i + 1u));
		_Matrix_Chain MC2(matrix_chain.begin() + (partition - i + 1u), matrix_chain.end());
		return MatrixChainMultiplication(MC1, S, i, partition) * MatrixChainMultiplication(MC2, S, partition + 1, j);
	}

	template<typename _Matrix_Chain>
	auto MatrixChainMultiplication(const _Matrix_Chain& matrix_chain) -> typename _Matrix_Chain::value_type
	{
		assert(matrix_chain.size() > 0u);
		size_t size = matrix_chain.size();
		matrix<size_s> S(size, size, 0);
		matrix<size_m> M(size, size, 0);
		std::pair<size_t, long long> partition{ 0,0 };
		for (size_t k = 0; k < size; ++k)
		{
			for (size_t i = 0; i + k < size; ++i)
			{
				partition = MatrixChainPartition(matrix_chain, M, i, i + k);
				S[i][i + k] = partition.first;
				M[i][i + k] = partition.second;
			}
		}
		std::cout << "dp-complexity:" << M[0][size - 1] << endl;
		return MatrixChainMultiplication(matrix_chain, S, 0, matrix_chain.size() - 1);
	}

	template<typename _Matrix_Chain>
	auto NormalMatrixChainMultiplication(const _Matrix_Chain& matrix_chain) -> typename _Matrix_Chain::value_type
	{
		if (matrix_chain.size() == 0) {
			return _Matrix_Chain::value_type();
		}
		if (matrix_chain.size() == 1) {
			return *matrix_chain.begin();
		}
		if (matrix_chain.size() == 2) {
			return *matrix_chain.begin() * *(matrix_chain.begin() + 1);
		}		
		_Matrix_Chain MC(matrix_chain.begin() + 1, matrix_chain.end());
		return *matrix_chain.begin() * NormalMatrixChainMultiplication(MC);
	}

} //namespace alg

#endif //_ALGORITHM_TEMPLATE_MATRIX_CHAIN_MULTIPLICATION_H_
