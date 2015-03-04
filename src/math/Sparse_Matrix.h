#pragma once
#include "sys.h"
#include <map>
struct Value_Description
{
	Value_Description() : value(0.0), row(0), column(0)
	{  }
	Value_Description(double v, uint32 r, uint32 c) : value(v), row(r), column(c)
	{  }
	double value;
	uint32 row;
	uint32 column;
};
// symmetric upper triangular row-oriented sparse matrix
//все индексы, которые пердаются в качестве параметров, начинаются с единицы
class Sparse_Matrix_SUR
{
public:
	Sparse_Matrix_SUR(uint32 nrows);
	~Sparse_Matrix_SUR();
	void start_training();
	void add_value(uint32 row, uint32 column, double value);
	void stop_training(bool copy_values_to_matrix=true);
	double& get_val_ij(uint32 row, uint32 column);
	double mult_vec_i(double *vec, uint32 n);
	bool search_index(uint32 row, uint32 column, uint32 &ind);
	void clear();
	void print_data();
	void print();
	void print_to_file(char *filename);
	void zero();
	void zero_block(uint32 n);
	double* get_values_array()
	{
		assert(n_values);
		return values;
	}
	uint32* get_columns_array()
	{
		assert(n_values);
		return columns;
	}
	uint32* get_iofeir_array()
	{
		assert(n_values);
		return iofeir;
	}
	uint32 get_n_values()
	{
		return n_values;
	}

	uint32 get_n_rows()
	{
		return n_rows;
	}
private:
	uint64 general_index(uint32 row, uint32 column);
	bool is_training;
	map<uint64, Value_Description> training_data;
	uint32 *n_in_row;
	uint32 n_rows; //количество строк и столбцов в квадратной матрице;...
	
	double dummy;
	//основные структуры данных для хранения разряженной матрицы (подробнее читай MKL manual)
	double *values;
	uint32 *columns;
	uint32 *iofeir; //index of first element in row
	uint32 n_values; //является параметром готовности работы с матрицей ( если > 0, то матрица готова к работе )
};


// row-oriented sparse matrix
//все индексы, которые пердаются в качестве параметров, начинаются с единицы
class Sparse_Matrix_R
{
public:
	Sparse_Matrix_R(uint32 nrows, uint32 ncolumns);
	~Sparse_Matrix_R();
	void start_training();
	void add_value(uint32 row, uint32 column, double value);
	void stop_training(bool copy_values_to_matrix=true);
	double& get_val_ij(uint32 row, uint32 column);
	double mult_vec_i(double *vec, uint32 n);
	double transpose_mult_vec_i(double *vec, uint32 n);
	void transpose_mult_vec(double *vec, double *res);
	bool search_index(uint32 row, uint32 column, uint32 &ind);
	void clear();
	void print_data();
	void print();
	void zero();
	double* get_values_array()
	{
		assert(n_values);
		return values;
	}
	uint32* get_columns_array()
	{
		assert(n_values);
		return columns;
	}
	uint32* get_iofeir_array()
	{
		assert(n_values);
		return iofeir;
	}
	uint32 get_n_values()
	{
		return n_values;
	}

	uint32 get_n_rows()
	{
		return n_rows;
	}
	uint32 get_n_columns()
	{
		return n_columns;
	}
private:
	uint64 general_index(uint32 row, uint32 column);
	bool is_training;
	map<uint64, Value_Description> training_data;
	uint32 *n_in_row;
	uint32 n_rows; //количество строк 
	uint32 n_columns; //и столбцов в прямоугольной матрице;...
	
	double dummy;
	//основные структуры данных для хранения разряженной матрицы (подробнее читай MKL manual)
	double *values;
	uint32 *columns;
	uint32 *iofeir; //index of first element in row
	uint32 n_values; //является параметром готовности работы с матрицей ( если > 0, то матрица готова к работе )
};
