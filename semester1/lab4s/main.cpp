#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>

#include "matr.h"
#include "pm.h"

int LoadFromFile( FILE *F, std::vector<std::vector<double>> &A, std::vector<double> &x0, double &Epsilon )
{
  int i, j, N;

  if (fscanf(F, "%i", &N) < 1)
    return 0;

  fscanf(F, "%lf", &Epsilon);

  A.clear();
  A.resize(N);
  x0.resize(N);

  for (i = 0; i < N; i++)
    A[i] = std::vector<double>(N);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      fscanf(F, "%lf", &(A[i][j]));

  for (i = 0; i < N; i++)
    fscanf(F, "%lf", &(x0[i]));

  return 1;
} /* End of 'LoadMatrixFromFile' function */

int main( void )
{
  std::vector<std::vector<double>> A;
  int i;
  std::vector<double> x0;
  eigen e;
  double eps;

  FILE *F = fopen("solve.out", "w");
  FILE *Fin = fopen("m.in", "r");

  if (F == NULL || Fin == NULL)
    return 0;

  while (LoadFromFile(Fin, A, x0, eps))
  {
    e = pmWithShift(A, x0, eps, 1.1);

    if (IsEqual(MatrMulVec(A, e.Vector), VecMulNum(e.Vector, e.Value), eps) == 0)
      printf("OK. Steps count: %i\n", e.Steps);
    else
      printf("BAD: %i\n", IsEqual(MatrMulVec(A, e.Vector), VecMulNum(e.Vector, e.Value), eps));

    fprintf(F, "%i\n", x0.size());
    fprintf(F, "%.16f\n", e.Value);
    for (i = 0; i < (int)e.Vector.size(); i++)
      fprintf(F, "%.16f ", e.Vector[i]);
    fprintf(F, "\n");
    fprintf(F, "%i\n\n", e.Steps);
  }

  fclose(F);
  fclose(Fin);

  return 0;
} /* End of 'main' function */
