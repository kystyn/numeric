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
void lieqsys::Diagonal( matr const &A, vec const &b, vec &x, uint N )
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
void lieqsys::Left( matr const &A, vec const &b, vec &x, uint N )
{
  x[0] = b[0] / A[0][0];

  for (uint k = 1; k < N; k++)
  {
    double p = 0;

    for (uint j = 0; j < k; j++)
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
void lieqsys::Right( matr const &A, vec const &b, vec &x, uint N )
{
  x[N - 1] = b[N - 1] / A[N - 1][N - 1];

  for (int k = N - 2; k >= 0; k--)
  {
    double p = 0;

    for (int j = k + 1; j < N; j++)
      p += A[k][j] * x[j];

    x[k] = (b[k] - p) / A[k][k];
  }
} /* End of 'Right' function */
