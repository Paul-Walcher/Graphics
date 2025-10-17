
#ifndef __MATRIX__
#define __MATRIX__

#include <vector>
#include <utility>
#include <iostream>
#include <functional>
#include <exception>
#include <cmath>
#include <tuple>


template<typename T>
class Matrix{

protected:

	int rows = 0;
	int columns = 0;

	std::vector<std::vector<T>> M;

	//0 element
	T zero_val;
	//1 element
	T one_val;

public:

	typedef T value_t;




	Matrix();
	Matrix(T zero_val, T one_val);
	Matrix(int rows, int columns, T init_value, T zero_val, T one_val);
	Matrix(int rows, int columns, std::vector<std::vector<T>> init, T zero_val, T one_val);

	static Matrix<T> identity(int size, T zero_val, T one_val);

	//shape of width and height
	const std::pair<int, int> shape() const;
	const std::pair<T, T> values() const;

	//get and set values
	T get(int row, int column) const;
	void set(int row, int column, T value);

	//get scalar value if 1x1 matrix
	T scalar() const;

	T sum() const;

	//inverse functions
	Matrix<T> inverse() const;
	Matrix<T> pseudo_inverse() const;

	//returns a rectangular slice of the matrix
	Matrix<T> slice(int rowstart, int rowend, int columnstart, int columnend) const;

	Matrix<T> transpose() const;

	//places the rectangular matrix inside this matrix with top left at (row, column)
	void place(int row, int column, const Matrix<T>& other);

	//for applying functions
	Matrix<T> apply(std::function<T(T)> f, int rowstart, int rowend, int columnstart, int columnend) const;

	Matrix<T> hstack(Matrix<T>& other);
	Matrix<T> vstack(Matrix<T>& other);

	void print() const;

	//uses the - operator on all elements
	Matrix<T> negate() const;

	Matrix<T> copy() const;

	T norm(int degree) const;

	Matrix<T> normalize();

	std::tuple<Matrix<T>, Matrix<T>> QR_decomposition() const;
	std::tuple<Matrix<T>, Matrix<T>> eigenpairs_QR() const;


	//operators
	Matrix<T> operator*(const T& value) const;
	Matrix<T> operator*(const Matrix<T>& other) const;
	Matrix<T>& operator*=(const T& value);
	Matrix<T>& operator*=(const Matrix<T>& other);

	Matrix<T> operator+(const T& value) const;
	Matrix<T> operator+(const Matrix<T>& other) const;
	Matrix<T>& operator+=(const T& value);
	Matrix<T>& operator+=(const Matrix<T>& other);

	bool operator==(const Matrix<T>& other) const;
	//cmp should be a function of value 1 and value 2, 
	//and return true if they are considered equal or acceptably close, false if not
	bool equals(const Matrix<T>& other, std::function<bool(T, T)> cmp) const;

	Matrix<T> operator[](std::pair<std::pair<int, int>, std::pair<int, int>> indices);


};


template<typename T>
Matrix<T>::Matrix(){
	
}

template<typename T>
Matrix<T>::Matrix(T zero_val, T one_val){
	this->zero_val = zero_val;
	this->one_val = one_val;
}

template<typename T>
Matrix<T>::Matrix(int rows, int columns, T init_value, T zero_val, T one_val){


	if (rows == 0 || columns == 0){
		throw std::invalid_argument("False construction of Matrix");
	}

	this->zero_val = zero_val;
	this->one_val = one_val;

	this->rows = rows;
	this->columns = columns;

	for(int i = 0; i < rows; i++){

		std::vector<T> V(columns, init_value);
		M.push_back(V);

	}

}

template<typename T>
Matrix<T>::Matrix(int rows, int columns, std::vector<std::vector<T>> init, T zero_val, T one_val){


	if (rows == 0 || columns == 0 || init.size() != rows || init[0].size() != columns){

		throw std::invalid_argument("False construction of Matrix");
	}

	this->zero_val = zero_val;
	this->one_val = one_val;

	this->rows = rows;
	this->columns = columns;

	M = init;

}

template<typename T>
Matrix<T> Matrix<T>::identity(int size, T zero_val, T one_val) {

	if (size < 0){
	
	throw std::invalid_argument("False construction of  Identity Matrix");
	
	}

	Matrix<T> mat(size, size, zero_val, zero_val, one_val);

	for(int i = 0; (i < size); i++){
		mat.set(i, i, one_val);
	}

	return mat;

}

template<typename T>
const std::pair<int, int> Matrix<T>::shape() const {

	std::pair<int, int> P;
	P.first = rows;
	P.second = columns;

	return P;

}

template<typename T>
T Matrix<T>::get(int row, int column) const {

	if (row < 0 || row >= rows || column < 0 || column >= columns){
		throw std::invalid_argument("Specified value does not exist in matrix");
	}

	return M[row][column];
}

template<typename T>
void Matrix<T>::set(int row, int column, T value){

	if (row < 0 || row >= rows || column < 0 || column >= columns){
		throw std::invalid_argument("Specified value does not exist in matrix");
	}

	M[row][column] = value;
}

template<typename T>
Matrix<T> Matrix<T>::slice(int rowstart, int rowend, int columnstart, int columnend) const {
	
	if (rowstart < 0 || rowend > rows || columnstart < 0 || columnend > columns 
			){
		throw std::invalid_argument("Arguments for slice are invalid");
	}

	Matrix<T> mat(rowend - rowstart, columnend - columnstart, zero_val, zero_val, one_val);

	for(int i = rowstart; i < rowend; i++){
		for(int j = columnstart; j < columnend; j++){
			mat.set(i - rowstart, j - columnstart, get(i, j));
		}
	}

	return mat;

}

template<typename T>
T Matrix<T>::scalar() const {

	if (!(columns == 1 && rows == 1)){
		throw std::runtime_error("Matrix is not in scalar form");
	}

	return get(0, 0);

}

template<typename T>
Matrix<T> Matrix<T>::copy() const {

	Matrix<T> back(rows, columns, zero_val, zero_val, one_val);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			back.set(i, j, get(i, j));

	return back;

}

template<typename T>
Matrix<T> Matrix<T>::hstack(Matrix<T>& other){

	if (other.rows != rows){
		throw std::invalid_argument("Matrices have to have the same number of rows.");
	}


	Matrix<T> back(rows, columns + other.columns, zero_val, zero_val, one_val);

	for (int i = 0; i < rows; i++){
		for(int j = 0; j < columns; j++){
			back.set(i, j, get(i, j));
		}
	}


	for (int i = 0; i < other.rows; i++){
		for(int j = 0; j < other.columns; j++){
			back.set(i, columns+j, other.get(i, j));
		}
	}



	return back;

}

template<typename T>
Matrix<T> Matrix<T>::vstack(Matrix<T>& other){

	if (other.columns != columns){
		std::cout << columns << " " << other.columns << "\n";
		throw std::invalid_argument("Matrices have to have the same number of columns");
	}

	Matrix<T> back(rows + other.rows, columns, zero_val, zero_val, one_val);

	for (int i = 0; i < rows; i++){
		for(int j = 0; j < columns; j++){
			back.set(i, j, get(i, j));
		}
	}


	for (int i = 0; i < other.rows; i++){
		for(int j = 0; j < other.columns; j++){
			back.set(rows+i, j, other.get(i, j));
		}
	}



	return back;

}

template<typename T>
void Matrix<T>::place(int row, int column, const Matrix<T>& other){

	if (row < 0 ||  column < 0 
		|| row + other.rows > rows || column + other.columns > columns
		){

		throw std::invalid_argument("Matrix does not fit at specified place");
	}

	for(int i = 0; i < other.rows; i++){
		for(int j = 0; j < other.columns; j++){
			set(row + i, column + j, other.get(i, j));
		}
	}

}


template<typename T>
void Matrix<T>::print() const {

	for(auto v : M){
		for ( auto k : v){
			std::cout << k << " ";
		}
		std::cout << "\n";
	}

}

template<typename T>
Matrix<T> Matrix<T>::negate() const {

	Matrix<T> other = (*this);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			other.set(i, j, -get(i, j));

	return other;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() const {


	Matrix<T> back(columns, rows, zero_val, zero_val, one_val);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			back.set(j, i, get(i, j));

	return back;
}

template<typename T>
T Matrix<T>::sum() const {

	T result = zero_val;

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < columns; j++){
			result += get(i, j);
		}
	}

	return result;

}

template<typename T>
T Matrix<T>::norm(int degree) const {

	if (!(rows == 1 || columns == 1)){
		throw std::invalid_argument("Matrix has to have vector form");
	}

	if (degree == -1){

		T v = zero_val;
		for(int i = 0; i < (std::max(rows, columns)); i++){

			int r = (rows == 0 ? 0 : i);
			int c = (columns == 1 ? 0 : i);

		
			v = std::max(v, get(r, c));
		}

		return v;

	}

	T sum = zero_val;

	for(int i = 0; i < (std::max(rows, columns)); i++){

		int r = (rows == 0 ? 0 : i);
		int c = (columns == 1 ? 0 : i);

		sum += std::abs(std::pow(get(r, c), degree));

	}

	return std::pow(sum, one_val / degree);

}

template<typename T>
Matrix<T> Matrix<T>::normalize(){

	return (*this) * (one_val / this->norm(2));

}

template<typename T>
Matrix<T> Matrix<T>::inverse() const {

	if (rows != columns){
		throw std::invalid_argument("Matrix has to be square!");
	}



	Matrix<T> I = Matrix<T>::identity(rows, zero_val, one_val);
	Matrix<T> A = (*this);


	//row reduction to stair form



	for(int i = 0; i < columns; i++){



		T svalue = A.get(i, i);

		std::function<T(T)> div = [&svalue](T v){return v * 1 / svalue;};

		A = A.apply(div, i, i+1, 0, columns);
		I = I.apply(div, i, i+1, 0, columns);

		//now subtracting
		Matrix<T> toprow = A[{{i, i+1}, {0, columns}}];
		Matrix<T> toprowi = I[{{i, i+1}, {0, columns}}];

		for(int j = i+1; j < rows; j++){


			Matrix<T> row = A[{{i,i+1}, {0, columns}}];

			if (row.sum() == 0.0){
				throw std::runtime_error("Matrix has no inverse");
			}

			svalue = A.get(j, i);

			Matrix<T> sub = (toprow * A.get(j, i)).negate();
			Matrix<T> subi = (toprowi * A.get(j, i)).negate();


			Matrix<T> crow = A[{{j, j+1}, {0, columns}}];

			Matrix<T> icrow = I[{{j, j+1}, {0, columns}}];


			crow += sub;
			icrow += subi;


			A.place(j, 0, crow);
			I.place(j, 0, icrow);

		}


	}


	//now eliminating the values behing the diagonal
	for(int i = rows-1; i >= 0; i--){


		Matrix<T> subrow = I[{{i, i+1}, {0, columns}}];

		for(int j = i-1; j >= 0; j--){

			T svalue = A.get(j, i);

			Matrix<T> crow = I[{{j, j+1}, {0, columns}}];

			crow += svalue * subrow * -1;

			I.place(j, 0, crow);

			A.set(j, i, zero_val);

		}

	}

	return I;

}

template<typename T>
Matrix<T> Matrix<T>::pseudo_inverse() const {

	Matrix<T> A = (*this);

	return (A.transpose() * A).inverse() * A.transpose();

}

//for applying functions
template<typename T>
Matrix<T> Matrix<T>::apply(std::function<T(T)> f, int rowstart, int rowend, int columnstart, int columnend) const {

	if (rowstart < 0 || rowend > rows || columnstart < 0 || columnend > columns 
			){
		throw std::invalid_argument("Arguments for slice are invalid");
	}

	Matrix<T> copy = (*this);

	for(int i = rowstart; i < rowend; i++){
		for(int j = columnstart; j < columnend; j++){

			copy.set(i, j, f(copy.get(i, j)));

		}
	}

	return copy;

}

template<typename T>
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::QR_decomposition() const {

	try{
		Matrix<T> I = inverse();
	}
	catch (std::exception& e){
		throw std::invalid_argument("The matrix has got to have a inverse!");
	}



	Matrix<T> Q(rows, columns, zero_val, zero_val, one_val);
	Matrix<T> R(rows, columns, zero_val, zero_val, one_val);

	std::vector<Matrix<T>> q;


	for(int i = 0; i < columns; i++){

		Matrix<T> a = slice(0, rows, i, i+1);
		Matrix<T> sum(rows, 1, zero_val, zero_val, one_val);


		for(int j = 0; j < i; j++){


			T v = (a.transpose() * q[j]).scalar();



			R.set(j, i, v);

			Matrix<T> nq = q[j] * v;

			sum += nq;

		}


		Matrix<T> qi = a + sum.negate();

		T normval = qi.norm(2);
		
		R.set(i, i, normval);

		qi = qi * (one_val / normval);

		Q.place(0, i, qi);

		q.push_back(qi);



	}



	return std::make_tuple(Q, R);

}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& value) const {


	Matrix<T> M2 = (*this);

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < columns; j++){
			M2.M[i][j] *= value;
		}
	}

	return (M2);

}

template<typename T>
Matrix<T> operator*(const typename Matrix<T>::value_t& value, const Matrix<T>& M) {

	return M * value;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {


	if (columns != other.rows){
		throw std::invalid_argument("M1 * M2: M1 must have the same number of columns as M2 has rows.");		
	} 



	Matrix<T> back(rows, other.columns, zero_val, zero_val, one_val);

	for(int i = 0; i < other.columns; i++){
		for(int j = 0; j < rows; j++){

			T result = zero_val;

			for(int k = 0; k < columns; k++){

				result += get(j, k) * other.get(k, i);

			}

			back.set(j, i, result);
		}
	}

	back.rows = rows;
	back.columns = other.columns;

	
	return back;

}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& value) {

	(*this) = (*this) * value;
	return (*this);

}


template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {



	(*this) = (*this) * other;
	return (*this);

}

//+

template<typename T>
Matrix<T> Matrix<T>::operator+(const T& value) const {


	Matrix<T> M2 = (*this);

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < columns; j++){
			M2.M[i][j] += value;
		}
	}

	return (M2);

}


template<typename T>
Matrix<T> operator+(const typename Matrix<T>::value_t& value, const Matrix<T>& M) {

	return M + value;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {


	if (columns != other.columns || rows != other.rows){
		throw std::invalid_argument("M1 * M2: M1 must have the same number of columns as M2 has rows.");		
	} 

	Matrix<T> back(rows, columns, zero_val, zero_val, one_val);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			back.set(i, j, get(i, j) + other.get(i, j));

	
	return back;

}



template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T& value){

	(*this) = (*this) + value;
	return (*this);

}


template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other){

	(*this) = (*this) + other;
	return (*this);

}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& other) const {


	if (columns != other.columns || rows != other.rows){
		return false;	
	} 

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			if (get(i, j) != other.get(i, j))
				return false;

	return true;

}

template<typename T>
bool Matrix<T>::equals(const Matrix<T>& other, std::function<bool(T, T)> cmp) const {


	if (columns != other.columns || rows != other.rows){
		return false;	
	} 

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			if (!(cmp(get(i, j), other.get(i, j))))
				return false;

	return true;

}

template<typename T>
Matrix<T> Matrix<T>::operator[](std::pair<std::pair<int, int>, std::pair<int, int>> indices){

	int x = indices.first.first;
	int w = indices.first.second;
	int y = indices.second.first;
	int h = indices.second.second;

	return slice(x, w, y, h);

}

template<typename T>

const std::pair<T, T> Matrix<T>::values() const{

	std::pair<T, T> p;
	p.first = zero_val;
	p.second = one_val;

	return p;

}

#endif