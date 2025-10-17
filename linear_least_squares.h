
#ifndef __LLS__
#define __LLS__

#include "matrix.h"

template<typename T>
Matrix<T> linear_least_squares(Matrix<T> A, Matrix<T> x, Matrix<T> b){

	std::pair<T, T> V = A.values();

	if (A.shape().first != b.shape().first){
		throw std::runtime_error("Shapes A and b are not aligned");
	}

	Matrix <T> x_back(A.shape().second, 1, V.first, V.first, V.second);

	for(int i = 0; i < x.shape().first; i++){

		T val = V.first;

		for(int j = 0; j < A.shape().first; j++){
			for(int k = 0; k < A.shape().second; k++){

				if (k == i) continue;

				val += x.get(k, 0) * A.get(j, k) * A.get(j, i);



			}
		}

		T div = V.first;

		for(int j = 0; j < A.shape().first; j++){
			div += std::pow(A.get(j, i), 2);
		}

		x_back.set(i, 0, val / div);

	}

	return x_back;

}

#endif