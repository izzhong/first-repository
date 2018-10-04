#include<iostream>
#include<vector>
#include<random>
#include<ctime>
#include<Windows.h>
#include"matrix_chain_multiplication.h"
using namespace std;
using namespace alg;

vector<matrix<double>> RandomMatrixChain(size_t max_row_col, size_t size)
{
	default_random_engine e(static_cast<size_t>(time(0)));
	uniform_int_distribution<size_t> ui(max_row_col);
	size_t row = ui(e);
	size_t col = ui(e);
	vector<matrix<double>> mc;
	for (size_t i = 0; i < size; ++i)
	{
		mc.push_back(matrix<double>(row, col,1.0));
		row = col;
		col = ui(e);
	}
	return mc;
}

int main()
{
	size_t max_row_col = 50;
	size_t size = 500;
	default_random_engine e(static_cast<size_t>(time(0)));
	uniform_int_distribution<size_t> ui(1,max_row_col);
	size_t row = ui(e);
	size_t col = ui(e);
	vector<matrix<double>> matrix_chain;
	for (size_t i = 0; i < size; ++i)
	{
		matrix_chain.emplace_back(row, col, 1.0);
		row = col;
		col = ui(e);
		std::cout << i << endl;
	}

	double frequency = 0.0;
	double start = 0.0;
	double end = 0.0;
	double normal_start = 0.0;
	double normal_end = 0.0;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);
	QueryPerformanceCounter((LARGE_INTEGER*)&normal_start);
	auto normal_result = NormalMatrixChainMultiplication(matrix_chain);
	QueryPerformanceCounter((LARGE_INTEGER*)&normal_end);

	QueryPerformanceCounter((LARGE_INTEGER*)&start);
	auto result = MatrixChainMultiplication(matrix_chain);
	QueryPerformanceCounter((LARGE_INTEGER*)&end);

	cout << "normal_complexity:" << MatrixMultiplicationComplexity(matrix_chain) << endl;
	cout << " numb :" << size << endl;
	cout <<"result:"<< (result == normal_result) << endl;
	cout <<"normal:"<< (normal_end - normal_start) / (frequency) << endl;
	cout << "  dp  :" << (end - start) / (frequency) << endl;
	cout << " rate :" << (normal_end - normal_start) / (end - start) << endl;

	system("pause");
	return 0;
}