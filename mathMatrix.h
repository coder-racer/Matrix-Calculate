#pragma once

class Matrix {
public:
	Matrix();

	int network[5][5];
	int GetStrings();
	int GetColumns();
	void SetStrings(const unsigned int strings);
	void SetColumns(const unsigned int columns);
	void Clear();
	static void AddPoinErrorFunc(void(*AddError)(const char* error));	
	static void multiplication(Matrix input, int num, Matrix* out);//умножение
	static void addition(Matrix input, Matrix input2, Matrix* out);//сложение
	static void subtraction(Matrix input, Matrix input2, Matrix* out);//вычитание
	static void matrixMultiplication(Matrix input, Matrix input2, Matrix* out);//перемножение
	static void minor(Matrix input, Matrix* out);//минор
	static void determinant(Matrix input, int* out);//детерминант 
	static void unMinor(Matrix input, Matrix* out);//алгоритмическое доп
	static void transpose(Matrix input, Matrix* out);
	static void inverseMatrix(Matrix input, Matrix* out, int * det);
private:
	static void determinant2(Matrix input, int* out);//детерминант - 2
	static void determinant3(Matrix input, int* out); //детерминант -3 
	static int GetShard(Matrix input, int row, int col);
	static int UnGetShard(Matrix input, int row, int col);
	static void (*AddError)(const char* error);
	unsigned int strings;
	unsigned int columns;
};
