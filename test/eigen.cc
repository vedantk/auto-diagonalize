#include <cmath>
#include <cstdio>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
	int fib = 42;

	MatrixXcd T;

	T = MatrixXcd(3, 3);

	MatrixXd F(2, 2);
	F(0, 0) = 0;
	F(0, 1) = 1;
	F(1, 0) = 1;
	F(1, 1) = 1;
	EigenSolver<MatrixXd> eigensolver(F);

	MatrixXcd D = eigensolver.eigenvalues().asDiagonal();
	D(0, 0) = pow(D(0, 0), fib);
	D(1, 1) = pow(D(1, 1), fib);

	MatrixXcd P = eigensolver.eigenvectors();
	MatrixXcd Pinv = P.inverse();
	MatrixXcd S = (P * D * Pinv);

	cout << F << endl
		<< P << endl
		<< Pinv << endl
		<< D << endl
		<< S << endl;

	cout << "fib(" << fib << "):\n";
	cout << S << endl;
	printf("%.0f\n", round(real(S(1, 1))));

	return 0;
}
