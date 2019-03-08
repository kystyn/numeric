#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "relax.h"
#include "matr.h"

int LoadFromFile( FILE *F, std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x0, double &Omega, double &Epsilon )
{
  int i, j, N;

  if (fscanf(F, "%i", &N) < 1)
    return 0;

  fscanf(F, "%lf%lf", &Omega, &Epsilon);

  A.clear();
  A.resize(N);
  b.resize(N);
  x0.resize(N);

  for (i = 0; i < N; i++)
    A[i] = std::vector<double>(N);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      fscanf(F, "%lf", &(A[i][j]));

  for (i = 0; i < N; i++)
    fscanf(F, "%lf", &(b[i]));

  for (i = 0; i < N; i++)
    fscanf(F, "%lf", &(x0[i]));

  return 1;
} /* End of 'LoadMatrixFromFile' function */

int main( void )
{
  std::vector<std::vector<double>> A;
  std::vector<double> b;
  int i, st;
  std::vector<double> x, x0;
  double omega, eps;

  FILE *F = fopen("solve.out", "w");
  FILE *Fin = fopen("m.inn", "r");

  if (F == NULL || Fin == NULL)
    return 0;

  while (LoadFromFile(Fin, A, b, x0, omega, eps))
  {
    x = Relax(A, b, x0, omega, eps, st);

    if (IsEqual(MatrMulVec(A, x), b, eps) == 0)
      printf("OK. Steps count: %i\n", st);
    else
      printf("BAD: %i\n", IsEqual(MatrMulVec(A, x), b, eps));

    fprintf(F, "%i\n", b.size());
    /*for (i = 0; i < b.size(); i++)
      fprintf(F, "%.16f ", b[i]);*/
    for (i = 0; i < b.size(); i++)
      fprintf(F, "%.16f ", x[i]);
    fprintf(F, "\n");
    fprintf(F, "%i\n", st);
  }

  fclose(F);
  fclose(Fin);

  return 0;
} /* End of 'main' function */
