#include "matrix.h"
#include "linear_least_squares.h"
#include <string>

typedef Matrix<float> fMatrix;

const float epsilon = 1e-5;

void assert(bool term, std::string msg){
	if (!term){
		std::cout << msg << "\n";
		exit(1);
	}
}

int main(){

	auto cmp = [](float a, float b){return std::abs(a - b) < epsilon;};

	//check every function

	fMatrix A(
				3, 3,
				{{1, 2, 3}, {5, 3, 1}, {9, 2, 3}},
				0, 1
		);

	fMatrix x(
				3, 1,
				{{5}, {5}, {3}},
				0, 1
		);

	fMatrix b(3, 1,
			{{8}, {3}, {4}},
			0, 1
		);




	std::cout << A.determinant() << "\n";


	

	return 0;
}



