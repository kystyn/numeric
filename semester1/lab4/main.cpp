#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>

#include "matr.h"
#include "pm.h"

bool LoadFromFile( FILE *F, matr &A, double &Epsilon, double &Shift )
{
  int i, j, N;

  if (fscanf(F, "%i", &N) < 1)
    return false;

  fscanf(F, "%lf", &Epsilon);
  fscanf(F, "%lf", &Shift);

  A = matr(N, N);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      fscanf(F, "%lf", &(A[i][j]));

  return true;
} /* End of 'LoadMatrixFromFile' function */

int main( void )
{
  matr A;
  vec x0;
  eigen e;
  double eps, shift;

  FILE *F = fopen("solve.out", "w");
  FILE *Fin = fopen("m.in", "r");

  if (F == NULL || Fin == NULL)
    return 0;

  while (LoadFromFile(Fin, A, eps, shift))
  {
    e = invItWithShift(A, eps, shift);

    fprintf(F, "%.16f\n", e.Value);
    fprintf(F, "%i\n\n", e.Steps);
  }

  fclose(Fin);
  fclose(F);

  return 0;
} /* End of 'main' function */
