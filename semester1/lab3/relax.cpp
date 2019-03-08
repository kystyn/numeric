#include "matr.h"

std::vector<double> Relax( std::vector<std::vector<double>> const &A, std::vector<double> const &b,
                           std::vector<double> const &x0, double Omega, double Epsilon, int &Steps ) {

  std::vector<double> t(A.size(), 0);
  std::vector<std::vector<double>> D, L = A, R = A;

  Steps = 0;

  for (int i = 0; i < A.size(); i++) {
    D.push_back(t);
    D[i][i] = A[i][i];
  }

  for (int y = 0; y < A.size(); y++)
    for (int x = 0; x < A.size(); x++) {
      if (x >= y)
        L[y][x] = 0;
      if (x <= y)
        R[y][x] = 0;
    }

  std::vector<std::vector<double>> x;

  x.push_back(x0);

  while (x.size() == 1 || (SqrNorm(VecAddVec(x.back(), VecMulNum(x[x.size() - 2], -1))) >= Epsilon * Epsilon)) {
    Steps++;
    for (int i = 0; i < b.size(); i++) {
      double s1 = 0, s2 = 0;

      for (int j = 0; j <= i - 1; j++)
        s1 += A[i][j] * t[j];

      for (int j = i + 1; j < b.size(); j++)
        s2 += A[i][j] * x.back()[j];

      t[i] = (1 - Omega) * x.back()[i] + Omega / A[i][i] * (b[i] - s1 - s2);
    }
    x.push_back(t);
  }

  return x.back();
} /* End of 'Relax' function */