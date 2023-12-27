/* Density Matrix Renormalization Group */
#include <Accelerate/Accelerate.h>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using Eigen ::MatrixXd;

/* Define Constant */
const int LatticeLength = 128;
const int MaxState = 32;
const double J = 1;
const double h = 1;
const double Sz[2][2] = {{1, 0}, {0, -1}};
const double Splus[2][2] = {{0, 1},
                            {0, 0}}; // S^- is just hermitian conj of S^+

int main(void) {
  std::cout << "Done!" << '\n';
  return 0;
}
