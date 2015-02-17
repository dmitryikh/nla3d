#include "Sparse_Matrix.h"

Sparse_Matrix_SUR::Sparse_Matrix_SUR(uint32 nrows) : n_rows(nrows)
{
	is_training = false;
	values = NULL;
	columns = NULL;
	n_in_row = NULL;
	n_values = 0;
	dummy = 0.0;
}

Sparse_Matrix_SUR::~Sparse_Matrix_SUR()
{
	clear();
}

void Sparse_Matrix_SUR::start_training()
{
	assert(n_rows);
	clear();
	if (is_training)
	{
		//удаляем данные о незаконченной тренировки и пишем warning
		warning("Sparse_Matrix_SUR::start_training: strart trainig before ending the previus one. The previus one has been deleted.");
	}
	is_training = true;
	if (n_in_row) delete[] n_in_row;
	n_in_row = new uint32[n_rows];
	memset(n_in_row,0,sizeof(uint32)*n_rows);
	training_data.clear();
	//сразу тренируем диагональные элементы! это нужно для решателя PARDISO
	for (uint32 i=0;i<n_rows;i++)
	{
		add_value(i+1,i+1,0.0);
	}
}

void Sparse_Matrix_SUR::add_value(uint32 row, uint32 column, double value)
{
	assert(row > 0 && row < n_rows+1);
	assert(column > 0 && column < n_rows+1);
	if (row > column) swap(column,row);
	uint32 gi = general_index(row, column);
	if (is_training)
	{
		if (training_data.count(gi) == 0)
		{
			training_data[gi] = Value_Description(value, row, column);
			n_in_row[row-1]++;
		}
		else
			training_data[gi].value += value;
	}
	else
	{
		uint32 index;
		if (search_index(row, column, index))
			values[index] += value;
		else
			warning("Sparse_Matrix_SUR::add_value: the position(%d,%d) is absent on current matrix", row,column);
	}
}

void Sparse_Matrix_SUR::stop_training(bool copy_values_to_matrix)
{
	n_values = training_data.size();
	values = new double[n_values];
	columns = new uint32[n_values];
	iofeir = new uint32[n_rows+1];
	iofeir[0] = 1;
	for (uint32 i=0; i < n_rows; i++)
	{
		iofeir[i+1] = iofeir[i]+n_in_row[i];
	}

	map<uint32, Value_Description>::iterator p = training_data.begin();
	uint32 counter = 0;
	while (p!=training_data.end())
	{
		//cout << p->first << "[" << p->second.row << "," << p->second.column << "]=  " << p->second.value << endl; 
		if (copy_values_to_matrix) 
			values[counter] = p->second.value;
		else
			values[counter] = 0.0;
		columns[counter] = p->second.column;
		p++;
		counter++;
	}
	assert(counter == n_values);
	//удаляем всю ненужную тренировочную информацию
	if (n_in_row) delete[] n_in_row, n_in_row = NULL;
	training_data.clear();
	is_training = false;
}

void Sparse_Matrix_SUR::clear()
{
	n_values = 0;
	if (values) delete[] values, values = NULL;
	if (columns) delete[] columns, columns = NULL;
	if (n_in_row) delete[] n_in_row, n_in_row = NULL;
}

void Sparse_Matrix_SUR::zero()
{
	assert(n_values);
	memset(values, 0, sizeof(double)*n_values);
}
uint32 Sparse_Matrix_SUR::general_index(uint32 row, uint32 column)
{
	assert(n_rows);
	//предпологается, что симметрия матрицы уже учтена (upper-triangular)
	//uint32 gen_index = 0;
	//uint32 k = 0;
	//while (k < row-1) gen_index += n_rows-k++;
	//gen_index += column - row;
	//return (uint32) //(row-2+n_rows)*(row-1)*0.5+column-row;
	double ind = ((double)n_rows+1.0)*((double)row-1.0) - (double)row*((double)row-1.0)*0.5+column-row;
	return (uint32) ind;
}
bool Sparse_Matrix_SUR::search_index(uint32 row, uint32 column, uint32 &ind)
{
	assert(n_values);
	assert(row > 0 && row < n_rows+1);
	assert(column > 0 && column < n_rows+1);
	if (row > column) swap(column,row);

	uint32 st = iofeir[row-1]-1;
	uint32 en = iofeir[row]-2;
	while(1)
	{
		if (en - st == 1)
		{
			if (columns[st] == column) 
			{
				ind = st;
				break;
			}
			if (columns[en] == column)
			{
				ind = en;
				break;
			}
			return false;
		}
		ind = (uint32) ((en+st)*0.5);
		
		if (columns[ind] == column) break;
		if (en == st)
			return false;

		if (columns[ind] > column) 
			en = ind;
		else
			st = ind;
	}
	return true;
}
double& Sparse_Matrix_SUR::get_val_ij(uint32 row, uint32 column)
{
	uint32 index;
	if (search_index(row, column, index))
		return values[index];
	else
		warning("Sparse_Matrix_SUR::get_val_ij: the position(%d,%d) is absent on current matrix", row,column);
	return dummy;
}


void Sparse_Matrix_SUR::print_data()
{
	cout << "values = {";
	for (uint32 i=0; i < n_values; i++)
		cout << values[i] << "   ";
	cout << "}" << endl;

	cout << "columns = {";
	for (uint32 i=0; i < n_rows; i++)
	{
		cout << "\n\n";
		for (uint32 j=iofeir[i]; j < iofeir[i+1]; j++)
			cout << columns[j-1] << " ";
	}
	cout << "}" << endl;

	cout << "iofeir = {";
	for (uint32 i=0; i < n_rows+1; i++)
		cout << iofeir[i] << "   ";
	cout << "}" << endl;
}

void Sparse_Matrix_SUR::print()
{
	uint32 ind;
	for (uint32 i=0; i < n_rows; i++)
	{
		cout << "[";
		for (uint32 j=0; j < n_rows;j++)
			if (search_index(i+1,j+1,ind))
				cout << values[ind] << "\t";
			else
				cout << "0\t";
		cout << "]" << endl;
	}
}
void Sparse_Matrix_SUR::print_to_file(char *filename)
{
	uint32 ind;
	ofstream file(filename);
	for (uint32 i=0; i < n_rows; i++)
	{
		for (uint32 j=0; j < n_rows;j++)
		{
			if (search_index(i+1,j+1,ind))
				file << values[ind];
			else
				file << "0";
			if (j != n_rows-1)
				file << "\t";
		}
		file << endl;
	}
	file.close();
}

double Sparse_Matrix_SUR::mult_vec_i(double *vec, uint32 n)
{
	double res = 0.0;
	double eps = 1e-20;
	uint32 index;


	for (uint32 j=0; j < n_rows; j++)
		if (fabs(vec[j]) > eps)
			if (search_index(n, j+1, index))
			{
				double ve = vec[j];
				double va = values[index];
				res += values[index]*vec[j];
			}
	return res;
}

void Sparse_Matrix_SUR::zero_block(uint32 n)
{
	assert(n > 0 && n < n_rows+1);
	assert(n_values);
	for (uint32 i=0; i < n; i++)
	{
		uint32 st = iofeir[i]-1;
		uint32 en = iofeir[i+1]-2;
		while (columns[en] > n)
			en--;
		memset(&values[st],0,sizeof(double)*(en-st+1)); //TODO: протестировать!
	}
}







//----------------------SPARSE_MATRIX_R----------------------------//

Sparse_Matrix_R::Sparse_Matrix_R(uint32 nrows, uint32 ncolumns) : n_rows(nrows), n_columns(ncolumns)
{
	is_training = false;
	values = NULL;
	columns = NULL;
	n_in_row = NULL;
	n_values = 0;
	dummy = 0.0;
}

Sparse_Matrix_R::~Sparse_Matrix_R()
{
	clear();
}

void Sparse_Matrix_R::start_training()
{
	assert(n_rows && n_columns);
	clear();
	if (is_training)
	{
		//удаляем данные о незаконченной тренировки и пишем warning
		warning("Sparse_Matrix_R::start_training: strart trainig before ending the previus one. The previus one has been deleted.");
	}
	is_training = true;
	n_in_row = new uint32[n_rows];
	memset(n_in_row,0,sizeof(uint32)*n_rows);
	training_data.clear();
}

void Sparse_Matrix_R::add_value(uint32 row, uint32 column, double value)
{
	assert(row > 0 && row < n_rows+1);
	assert(column > 0 && column < n_columns+1);
	uint32 gi = general_index(row, column);
	if (is_training)
	{
		if (training_data.count(gi) == 0)
		{
			training_data[gi] = Value_Description(value, row, column);
			n_in_row[row-1]++;
		}
		else
			training_data[gi].value += value;
	}
	else
	{
		uint32 index;
		if (search_index(row, column, index))
			values[index] += value;
		else
			warning("Sparse_Matrix_R::add_value: the position(%d,%d) is absent on current matrix", row,column);
	}
}

void Sparse_Matrix_R::stop_training(bool copy_values_to_matrix)
{
	n_values = training_data.size();
	values = new double[n_values];
	columns = new uint32[n_values];
	iofeir = new uint32[n_rows+1];
	iofeir[0] = 1;
	for (uint32 i=0; i < n_rows; i++)
	{
		iofeir[i+1] = iofeir[i]+n_in_row[i];
	}
	//iofeir[n_rows] = iofeir[n_rows-1] + 1;

	map<uint32, Value_Description>::iterator p = training_data.begin();
	uint32 counter = 0;
	while (p!=training_data.end())
	{
		if (copy_values_to_matrix) 
			values[counter] = p->second.value;
		else
			values[counter] = 0.0;
		columns[counter] = p->second.column;
		p++;
		counter++;
	}
	assert(counter == n_values);
	//удаляем всю ненужную тренировочную информацию
	if (n_in_row) delete[] n_in_row, n_in_row = NULL;
	training_data.clear();
	is_training = false;
}

void Sparse_Matrix_R::clear()
{
	n_values = 0;
	if (values) delete[] values, values = NULL;
	if (columns) delete[] columns, columns = NULL;
	if (n_in_row) delete[] n_in_row, n_in_row = NULL;
}

void Sparse_Matrix_R::zero()
{
	assert(n_values);
	memset(values, 0, sizeof(double)*n_values);
}
uint32 Sparse_Matrix_R::general_index(uint32 row, uint32 column)
{
	assert(n_columns);
	return (uint32) (row-1)*n_columns+column - 1;
}
bool Sparse_Matrix_R::search_index(uint32 row, uint32 column, uint32 &ind)
{
	assert(n_values);
	assert(row > 0 && row < n_rows+1);
	assert(column > 0 && column < n_columns+1);

	if (iofeir[row]-iofeir[row-1] == 0)
		return false;

	uint32 st = iofeir[row-1]-1;
	uint32 en = iofeir[row]-2;
	while(1)
	{
		if (en - st == 1)
		{
			if (columns[st] == column) 
			{
				ind = st;
				break;
			}
			if (columns[en] == column)
			{
				ind = en;
				break;
			}
			return false;
		}
		ind = (uint32) ((en+st)*0.5);
		
		if (columns[ind] == column) break;
		if (en == st)
			return false;

		if (columns[ind] > column) 
			en = ind;
		else
			st = ind;
	}
	return true;
}
double& Sparse_Matrix_R::get_val_ij(uint32 row, uint32 column)
{
	uint32 index;
	if (search_index(row, column, index))
		return values[index];
	else
		warning("Sparse_Matrix_R::get_val_ij: the position(%d,%d) is absent on current matrix", row,column);
	return dummy;
}


void Sparse_Matrix_R::print_data()
{
	cout << "values = {";
	for (uint32 i=0; i < n_values; i++)
		cout << values[i] << "   ";
	cout << "}" << endl;

	cout << "columns = {";
	for (uint32 i=0; i < n_rows; i++)
	{
		cout << "\n\n";
		for (uint32 j=iofeir[i]; j < iofeir[i+1]; j++)
			cout << columns[j-1] << " ";
	}
	cout << "}" << endl;

	cout << "iofeir = {";
	for (uint32 i=0; i < n_rows+1; i++)
		cout << iofeir[i] << "   ";
	cout << "}" << endl;
}

void Sparse_Matrix_R::print()
{
	uint32 ind;
	//ofstream file("mat.txt");
	for (uint32 i=0; i < n_rows; i++)
	{
		cout << "[";
		for (uint32 j=0; j < n_columns;j++)
			if (search_index(i+1,j+1,ind))
				cout << values[ind] << "\t";
			else
				cout << "0\t";
		cout << "]" << endl;
	}
	//file.close();
}

double Sparse_Matrix_R::mult_vec_i(double *vec, uint32 n)
{
	double res = 0.0;
	double eps = 1e-20;

	if (iofeir[n]-iofeir[n-1] == 0)
		return 0.0; //нет элементов в данной строке
	uint32 st = iofeir[n-1]-1;
	uint32 en = iofeir[n]-2;

	for (uint32 j=st; j <= en; j++)
				res += values[j]*vec[columns[j]-1];
	return res;
}

double Sparse_Matrix_R::transpose_mult_vec_i(double *vec, uint32 n)
{
	double res = 0.0;
	double eps = 1e-20;
	uint32 index;
	for (uint32 j=0; j < n_rows; j++)
		if (fabs(vec[j]) > eps)
			if (search_index(j+1, n, index))
			{
				double ve = vec[j];
				double va = values[index];
				res += values[index]*vec[j];
			}
	return res;
}

void Sparse_Matrix_R::transpose_mult_vec(double *vec, double *res)
{
	assert(vec && res);
	memset(res,0,sizeof(double)*n_columns);
	double eps = 1e-20;

	for (uint32 i=0; i < n_rows; i++)
	{
		if (iofeir[i+1]-iofeir[i] == 0)
			continue;
		uint32 st = iofeir[i]-1;
		uint32 en = iofeir[i+1]-2;
		for (uint32 j=st; j <= en; j++)
				res[columns[j]-1] += values[j]*vec[i];
	}
}
