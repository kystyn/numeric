#include <stdlib.h>
#include <math.h>
#include "matr.h"

/* Matrix mul vector function.
 * ARGUMENTS:
 *   matrix (M x N);
 *   vector (N);
 *   M, N;
 * RETURNS: vector (M).
 */
double * MatrMulVec( double **A, double *b, int M, int N )
{
  int i, j;
  double *x = malloc(sizeof(double) * N);

  for (i = 0; i < M; i++)
  {
    double r = 0;
    for (j = 0; j < N; j++)
      r += A[i][j] * b[j];
    x[i] = r;
  }

  return x;
} /* End of 'MatrMulVec' function */

/* Vector mul matrix function.
 * ARGUMENTS:
 *   vector (transposed) (M);
 *   matrix (M x N);
 *   M, N;
 * RETURNS: vector (transposed) (N).
 */
double * VecMulMatr( double *b, double **A, int M, int N )
{
  int i, j;
  double *x = malloc(sizeof(double) * N);

  for (i = 0; i < M; i++)
  {
    double r = 0;
    for (j = 0; j < N; j++)
      r += b[i] * A[j][i];
    x[i] = r;
  }

  return x;
} /* End of 'VecMulMatr' function */

 /* Matrix mul matrix function.
 * ARGUMENTS:
 *   matrix (M x N), (N X K);
 *   M, N, K;
 * RETURNS: matrix (M X K).
 */
double ** MatrMulMatr( double **A, double **B, int M, int N, int K )
{
  double **C = malloc(sizeof(double *) * M);
  int i, j, k;

  for (i = 0; i < M; i++)
    C[i] = malloc(sizeof(double) * K);

  for (i = 0; i < M; i++)
    for (j = 0; j < K; j++)
    {
      C[i][j] = 0;
      for (k = 0; k < N; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return C;
} /* End of 'MatrMulMatr' function */

/* Transpose matrix function.
 * ARGUMENTS:
 *   - matrix M X N;  
 *   - M, N;
 * RETURNS: matrix N X M.
 */
double ** Transposing( double **A, int M, int N )
{
  double **At = malloc(sizeof(double *) * N);
  int i, j;

  for (i = 0; i < N; i++)
    At[i] = malloc(sizeof(double) * M);

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      At[i][j] = A[j][i];

  return At;
} /* End of 'Transposing' function */

/* Vec mul num function.
 */
double * VecMulNum( double *x, int N, double Num )
{
  double *v = malloc(sizeof(double) * N);
  int i;

  for (i = 0; i < N; i++)
    v[i] = x[i] * Num;

  return v;
} /* End of 'VecMulNum' function */

/* Dot product function.
 * ARGUMENTS:
 *   - vectors:
 *       X, Y (N);
 * RETURNS: dot product.
 */
double DotProduct( double *X, double *Y, int N )
{
  int i;
  double dp = 0;

  for (i = 0; i < N; i++)
    dp += X[i] * Y[i];

  return dp;
} /* End of 'DotProduct' function */

/* Compare two vectors function.
 * ARGUMENTS:
 *   vectors X(N), Y(N);
 * RETURNS: 1 if equal.
 */
int IsEqual( double *X, double *Y, int N )
{
  int i;
  double eps = 1e-16;
  double *a = malloc(sizeof(double) * N);

  for (i = 0; i < N; i++)
    a[i] = X[i] - Y[i];

  if (DotProduct(a, a, N) < eps)
    return 0;
  return (int)(log(DotProduct(a, a, N)) / log(10.0) + 0.5);
} /* End of 'IsEqual' function */
