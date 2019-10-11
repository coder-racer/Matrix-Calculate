#include "mathMatrix.h"

void Matrix::multiplication(Matrix input, int num, Matrix* out)
{
	out->SetColumns(input.GetColumns());
	out->SetStrings(input.GetStrings());
	for (int i = 0; i < input.GetStrings(); i++) {
		for (int j = 0; j < input.GetColumns(); j++) {
			out->network[i][j] = input.network[i][j] * num;
		}
	}
}

void Matrix::addition(Matrix input, Matrix input2, Matrix* out)
{
	if (input.GetColumns() != input2.GetColumns() || input.GetStrings() != input2.GetStrings())
	{
		Matrix::AddError(u8"Размерность матриц не совпадает.");
		return;
	}

	out->SetColumns(input2.GetColumns());
	out->SetStrings(input2.GetStrings());
	for (int i = 0; i < input.GetStrings(); i++) {
		for (int j = 0; j < input.GetColumns(); j++) {
			out->network[i][j] = input.network[i][j] + input2.network[i][j];
		}
	}
}

void Matrix::subtraction(Matrix input, Matrix input2, Matrix* out)
{
	if (input.GetColumns() != input2.GetColumns() || input.GetStrings() != input2.GetStrings())
	{
		Matrix::AddError(u8"Размерность матриц не совпадает.");
		return;
	}
	out->SetColumns(input2.GetColumns());
	out->SetStrings(input2.GetStrings());
	for (int i = 0; i < input.GetStrings(); i++) {
		for (int j = 0; j < input.GetColumns(); j++) {
			out->network[i][j] = input.network[i][j] - input2.network[i][j];
		}
	}
}

void Matrix::matrixMultiplication(Matrix input, Matrix input2, Matrix* out)
{
	if (input.GetColumns() != input2.GetStrings())
	{
		Matrix::AddError(u8"Столбцы первой матрицы не равны строкам второй.");
		return;
	}
	int size = input.GetColumns();
	out->SetColumns(size);
	out->SetStrings(size);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			out->network[i][j] = 0;

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int repeat = 0; repeat < size; repeat++)
			{
				out->network[i][j] += input.network[i][repeat] * input2.network[repeat][j];
			}
		}
	}
}

int Matrix::GetShard(Matrix input, int row, int col)
{
	Matrix temp;
	int size = input.GetColumns() - 1;
	temp.SetColumns(size);
	temp.SetStrings(size);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
			temp.network[i][j] = input.network[i][j];
		for (int j = col; j < size; j++)
			temp.network[i][j] = input.network[i][j + 1];
	}

	for (int i = row; i < size; i++)
	{
		for (int j = col; j < size; j++)
			temp.network[i][j] = input.network[i + 1][j + 1];
		for (int j = 0; j < col; j++)
			temp.network[i][j] = input.network[i + 1][j];
	}
	int il = 0;
	Matrix::determinant(temp, &il);
	return il;
}

int Matrix::UnGetShard(Matrix input, int row, int col)
{
	int index = row + col;
	int Shard = Matrix::GetShard(input, row, col);
	if (index % 2)
		Shard = -1 * Shard;
	return Shard;
}

void Matrix::minor(Matrix input, Matrix* out)
{
	if (input.GetColumns() != input.GetStrings())
	{
		Matrix::AddError(u8"Матрица должна быть квадратной.");
		return;
	}
	int size = input.GetStrings();
	out->SetColumns(size);
	out->SetStrings(size);
	for (int i = 0; i < size; i++)
	{
		for (int j= 0; j < size; j++)
		{
			out->network[i][j] = Matrix::GetShard(input, i, j);
		}
	}

}

void Matrix::determinant(Matrix input, int* out)
{
	if (input.GetColumns() != input.GetStrings())
	{
		Matrix::AddError(u8"Матрица должна быть квадратной.");
		return;
	}
	int size = input.GetStrings();
	if (size == 2)
	{
		Matrix::determinant2(input, out);
		return;
	}
	else if (size == 3)
	{
		Matrix::determinant3(input, out);
		return;
	}
	else
	{
		int F = 0;
		for (int i = 0; i < size; i++)
		{
			F += input.network[0][i] * Matrix::UnGetShard(input, 0, i);
		}
		*out = F;
	}

}

void Matrix::unMinor(Matrix input, Matrix* out)
{
	if (input.GetColumns() != input.GetStrings())
	{
		Matrix::AddError(u8"Матрица должна быть квадратной.");
		return;
	}
	int size = input.GetStrings();
	out->SetColumns(size);
	out->SetStrings(size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out->network[i][j] = Matrix::UnGetShard(input, i, j);
		}
	}
}

void Matrix::transpose(Matrix input, Matrix* out)
{
	if (input.GetColumns() != input.GetStrings())
	{
		Matrix::AddError(u8"Матрица должна быть квадратной.");
		return;
	}
	int size = input.GetStrings();
	out->SetColumns(size);
	out->SetStrings(size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out->network[i][j] = input.network[j][i];
		}
	}

}

void Matrix::inverseMatrix(Matrix input, Matrix* out, int* det)
{
	*det = 0;
	Matrix::determinant(input, det);
	Matrix trans;
	Matrix::transpose(input, &trans);
	Matrix::unMinor(trans, out);
}

void Matrix::determinant2(Matrix input, int* out)
{
	*out = input.network[0][0] * input.network[1][1] - input.network[0][1] * input.network[1][0];
}

/*
	[0][0] [0][1] [0][2]|[0][0] [0][1]
						|
	[1][0] [1][1] [1][2]|[1][0] [1][1]
						|
	[2][0] [2][1] [2][2]|[2][0] [2][1]
*/
void Matrix::determinant3(Matrix input, int* out)
{	
	*out = input.network[0][0] * input.network[1][1] * input.network[2][2]
		+ input.network[0][1] * input.network[1][2] * input.network[2][0]
		+ input.network[0][2] * input.network[1][0] * input.network[2][1]
		- input.network[2][0] * input.network[1][1] * input.network[0][2]
		- input.network[2][1] * input.network[1][2] * input.network[0][0]
		- input.network[2][2] * input.network[1][0] * input.network[0][1]
		;

}

void Matrix::Clear()
{
	this->strings = 2;
	this->columns = 2;
	for (int row = 0; row < 5; row++)
	{
		for (int col = 0; col < 5; col++)
		{
			this->network[row][col] = 0;
		}
	}
}

int Matrix::GetStrings()
{
	return strings;
};

int Matrix::GetColumns()
{
	return columns;
};

void Matrix::SetStrings(const unsigned int strings)
{
	this->strings = strings;
}

void Matrix::SetColumns(const unsigned int columns)
{
	this->columns = columns;
}

Matrix::Matrix()
{
	this->Clear();
}

void Matrix::AddPoinErrorFunc(void(*AddError)(const char* error))
{
	Matrix::AddError = AddError;
}

void (*Matrix::AddError)(const char* error) = nullptr;
