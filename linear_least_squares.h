
#ifndef __LLS__
#define __LLS__

#include "matrix.h"
#include <iostream>

//loss function

typedef Matrix<float> fMatrix;

float loss(fMatrix X, fMatrix Y){

	fMatrix C = X + Y.negate();
	C = C.hadamard(C);

	return C.sum();

}

fMatrix dY(fMatrix X, fMatrix Y, float N){

	return 2 * (X + Y.negate()) * (1 / N);

}



void test(){

int m = 2;
int n = 5;
int p = 4;

fMatrix W = fMatrix::frand(m, n, 0.1, 1.0);
fMatrix B = fMatrix::frand(m, 1, 0.1, 1.0);


fMatrix X(4, 5, 
{

	{4, 4, 2.00, 0, 0},
	{3.8, 3.9, 2.55, 0, 0},
	{1.0, 5.0, 0, 2.55, 0},
	{2.0, 6.0, 0, 1.80, 0}

},
0.0,
1.0
	);

X = X.transpose();


fMatrix Y(2, 4, 

	{{1.0, 1.0, 0.0, 0.0},
	{0.0, 0.0, 1.0, 1.0}}
,
0.0, 
1.0

	);


	float lr = 0.1;


	for(int _ = 0; _ < 100; _++){

		fMatrix res = W * X;



		for(int i = 0; i < m; i++){
			for (int j = 0; j < p; j++){


				res.set(i, j, res.get(i, j) + B.get(i, 0));
			}
		}


		std::cout << loss(res, Y) << "\n";

		//gradient descent

		fMatrix dres = dY(res, Y, n*p);
		W += (-lr) * (dres * X.transpose());

		for(int i = 0; i < m; i++){
			float sum = 0.0;
			for(int j = 0; j < p; j++){
				sum += dres.get(i, j);
			}

			B.set(i, 0, B.get(i, 0) - (lr * sum));
		}


	}


	std::cout << "E" << "\n";

	fMatrix x(1, 5, {{2, 2, 1.50, 0.50, 0.0}}, 0, 1);

	x = x.transpose();

	fMatrix res = W * x;

	p = 1;


	for(int i = 0; i < m; i++){
		for (int j = 0; j < p; j++){


			res.set(i, j, res.get(i, j) + B.get(i, 0));
		}
	} 

	res.print();

}



#endif