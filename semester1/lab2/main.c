#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ldr.h"
#include "matr.h"

int LoadFromFile( FILE *F, double ***A, double **b, int *N )
{
  int i, j;

  if (fscanf(F, "%i", N) < 1)
    return 0;

  *A = malloc(*N * sizeof(double *));

  for (i = 0; i < *N; i++)
    (*A)[i] = malloc(sizeof(double) * *N);

  for (i = 0; i < *N; i++)
    for (j = 0; j < *N; j++)
      fscanf(F, "%lf", &((*A)[i][j]));

  *b = malloc(sizeof(double) * *N);

  for (i = 0; i < *N; i++)
    fscanf(F, "%lf", &((*b)[i]));

  return 1;
} /* End of 'LoadMatrixFromFile' function */

int main( void )
{
  double **A = NULL, *b = NULL;
  int N, i;
  double *x;

  FILE *F = fopen("solve.out", "w");
  FILE *Fin = fopen("m.in", "r");

  if (F == NULL || Fin == NULL)
    return;

  while (LoadFromFile(Fin, &A, &b, &N))
  {
    x = malloc(sizeof(double) * N);
    LDRSolve(A, b, x, N);

    if (IsEqual(MatrMulVec(A, x, N, N), b, N) == 0)
      printf("OK\n");
    else
      printf("BAD: %i\n", IsEqual(MatrMulVec(A, x, N, N), b, N));

    for (i = 0; i < N; i++)
      fprintf(F, "%.16f ", b[i]);
    fprintf(F, "\n");
    for (i = 0; i < N; i++)
      fprintf(F, "%.16f ", x[i]);
    fprintf(F, "\n\n");
  }

  fclose(F);
  fclose(Fin);

  return 0;
} /* End of 'main' function */
