#include <stdlib.h>

#include "ldr.h"
#include "simple.h"
#include "matr.h"

/* LDR decomposition function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - non-allocated L, D, R matrices;
 * RETURNS: None.
 */
void LDRDecomposition( double **A, double ***L, double ***D, double ***R, int N )
{
  int i, k;

  /* Variables for LDR */
  double
    *f = malloc(sizeof(double) * N),
    *g = malloc(sizeof(double) * N),
    *x = malloc(sizeof(double) * N),
    *y = malloc(sizeof(double) * N);
  double beta;

  *L = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*L)[i] = malloc(sizeof(double) * N);

  *D = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*D)[i] = malloc(sizeof(double) * N);

  *R = malloc(sizeof(double *) * N);
  for (i = 0; i < N; i++)
    (*R)[i] = malloc(sizeof(double) * N);

  /* A1 = L1 * D1 * R1 */

  (*L)[0][0] = 1;
  (*D)[0][0] = A[0][0];
  (*R)[0][0] = 1;

  for (k = 2; k <= N; k++)
  {
    /* Fill g_k-1T */
    for (i = 0; i < k - 1; i++)
      g[i] = A[k - 1][i];
    /* Fill f_k-1 */
    for (i = 0; i < k - 1; i++)
      f[i] = A[i][k - 1];

    /*** Set Lk ***/
    /* Solve X */
    Left(MatrMulMatr(Transposing(*R, k - 1, k - 1), Transposing(*D, k - 1, k - 1), k - 1, k - 1, k - 1),
      g, x, k - 1);

    for (i = 0; i < k - 1; i++)
    {
      (*L)[k - 1][i] = x[i];
      (*L)[i][k - 1] = 0;
    }
    (*L)[k - 1][k - 1] = 1;

    /*** Set Rk ***/
    /* Solve Y */
    Left(MatrMulMatr(*L, *D, k - 1, k - 1, k - 1), f, y, k - 1);

    for (i = 0; i < k - 1; i++)
    {
      (*R)[i][k - 1] = y[i];
      (*R)[k - 1][i] = 0;
    }
    (*R)[k - 1][k - 1] = 1;

    /*** Set Dk ***/
    /* Solve beta */
    beta = A[k - 1][k - 1] - DotProduct(VecMulMatr(x, *D, k - 1, k - 1), y, k - 1);
    for (i = 0; i < k - 1; i++)
    {
      (*D)[i][k - 1] = 0;
      (*D)[k - 1][i] = 0;
    }
    (*D)[k - 1][k - 1] = beta;
  }

  free(f);
  free(g);
  free(x);
  free(y);
} /* End of 'LDRDecomposition' function */

/* LDR solve function.
 * ARGUMENTS:
 *   - matrix (N X N);
 *   - right vector b (N);
 *   - N;
 *   - allocated result vector x (N);
 * RETURNS: None.
 */
void LDRSolve( double **A, double *b, double *x, int N )
{
  double **L, **D, **R;
  double
    *y = malloc(sizeof(double) * N),
    *z = malloc(sizeof(double) * N);

  LDRDecomposition(A, &L, &D, &R, N);

  /* Ly = b */
  Left(L, b, y, N);
  /* Dz = y */
  Diagonal(D, y, z, N);
  /* Rx = z */
  Right(R, z, x, N);

  free(y);
  free(z);
} /* End of 'LDRSolve' function */