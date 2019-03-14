#include <stdlib.h>
#include "simple.h"

using namespace mth;

/* Diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> &A;
 *   - size:
 *       int N;
 *   - right vector:
 *       std::vector<double> &b;
 *   - result (allocated) vector:
 *       std::vector<double> &x;
 * RETURNS: None.
 */
void lieqsys::Diagonal( matr const &A, vec const &b, vec &x, int N )
{
  for (int k = 0; k < N; k++)
    x[k] = b[k] / A[k][k];
} /* End of 'Diagonal' function */

/* Left diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> &A;
 *   - size:
 *       int N;
 *   - right vector:
 *       std::vector<double> &b;
 *   - result (allocated) vector:
 *       std::vector<double> &x;
 * RETURNS: None.
 */
void lieqsys::Left( matr const &A, vec const &b, vec &x, int N )
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
 *       std::vector<std::vector<double>> &A;
 *   - size:
 *       int N;
 *   - right vector:
 *       std::vector<double> &b;
 *   - result (allocated) vector:
 *       std::vector<double> &x;
 * RETURNS: None.
 */
void lieqsys::Right( matr const &A, vec const &b, vec &x, int N )
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
