#include <stdlib.h>
#include "simple.h"

/* Diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Diagonal( double *A[], double *b, double *x, int N )
{
  int k;

  for (k = 0; k < N; k++)
    x[k] = b[k] / A[k][k];

} /* End of 'Diagonal' function */

/* Left diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Left( double **A, double *b, double *x, int N )
{
  int k;

  x[0] = b[0] / A[0][0];

  for (k = 1; k < N; k++)
  {
    int j;
    double p = 0;

    for (j = 0; j < k; j++)
      p += A[k][j] * x[j];

    x[k] = (b[k] - p) / A[k][k];
  }
} /* End of 'Left' function */

/* Right diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Right( double *A[], double *b, double *x, int N )
{
  int k;

  x[N - 1] = b[N - 1] / A[N - 1][N - 1];

  for (k = N - 2; k >= 0; k--)
  {
    int j;
    double p = 0;

    for (j = k + 1; j < N; j++)
      p += A[k][j] * x[j];

    x[k] = (b[k] - p) / A[k][k];
  }
} /* End of 'Right' function */