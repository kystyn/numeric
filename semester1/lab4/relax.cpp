#include "matr.h"

vec Relax( matr const &A, vec const &b, vec const &x0, double Omega, double Epsilon, int &Steps ) {

  vec t(A.getH(), 0);

  matr D = A, L = A, R = A;

  Steps = 0;

  for (int i = 0; i < A.getH(); i++) {
    D[i] = std::vector<double>(A.getH(), 0);
    D[i][i] = A[i][i];
  }

  for (int y = 0; y < A.getH(); y++)
    for (int x = 0; x < A.getH(); x++) {
      if (x >= y)
        L[y][x] = 0;
      if (x <= y)
        R[y][x] = 0;
    }

  vec x = x0, prevx;

  while (x.cmp(x0, Epsilon) == matr::OK || !(x - prevx) >= Epsilon) {
    Steps++;
    for (int i = 0; i < b.getN(); i++) {
      double s1 = 0, s2 = 0;

      for (int j = 0; j <= i - 1; j++)
        s1 += A[i][j] * t[j];

      for (int j = i + 1; j < b.getN(); j++)
        s2 += A[i][j] * x[j];

      t[i] = (1 - Omega) * x[i] + Omega / A[i][i] * (b[i] - s1 - s2);
    }
    prevx = x;
    x = t;
  }

  return x;
} /* End of 'Relax' function */