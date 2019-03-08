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
std::vector<double> MatrMulVec( std::vector<std::vector<double>> &A, std::vector<double> &b )
{
  int i, j;
  std::vector<double> x(A[0].size());

  for (i = 0; i < b.size(); i++)
  {
    double r = 0;
    for (j = 0; j < A[0].size(); j++)
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
std::vector<double> VecMulMatr( std::vector<double> &b, std::vector<std::vector<double>> &A )
{
  int i, j;
  std::vector<double> x = std::vector<double>(A[0].size());

  for (i = 0; i < b.size(); i++)
  {
    double r = 0;
    for (j = 0; j < A[0].size(); j++)
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
std::vector<std::vector<double>> MatrMulMatr( std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B )
{
  std::vector<std::vector<double>> C = std::vector<std::vector<double>>(A.size());
  int i, j, k;

  for (i = 0; i < A.size(); i++)
    A[i].resize(B[0].size());

  for (i = 0; i < A.size(); i++)
    for (j = 0; j < B[0].size(); j++)
    {
      C[i][j] = 0;
      for (k = 0; k < A[0].size(); k++)
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
std::vector<std::vector<double>> Transposing( std::vector<std::vector<double>> &A )
{
  std::vector<std::vector<double>> At = std::vector<std::vector<double>>(A[0].size());
  int i, j;

  for (i = 0; i < At.size(); i++)
    At[i].resize(A.size());

  for (i = 0; i < A[0].size(); i++)
    for (j = 0; j < A.size(); j++)
      At[i][j] = A[j][i];

  return At;
} /* End of 'Transposing' function */

/* Vec mul num function.
 */
std::vector<double> VecMulNum( std::vector<double> &x, double Num )
{
  std::vector<double> v = std::vector<double>(x.size());
  int i;

  for (i = 0; i < x.size(); i++)
    v[i] = x[i] * Num;

  return v;
} /* End of 'VecMulNum' function */

/* Vec add vec  function.
 */
std::vector<double> VecAddVec( std::vector<double> &x, std::vector<double> &y ) {

  std::vector<double> p;

  for (int i = 0; i < x.size(); i++)
    p.push_back(x[i] + y[i]);

  return p;
} /* End of 'VecAddVec' function */

/* Dot product function.
 * ARGUMENTS:
 *   - vectors:
 *       X, Y (N);
 * RETURNS: dot product.
 */
double DotProduct( std::vector<double> &X, std::vector<double> &Y )
{
  int i;
  double dp = 0;

  for (i = 0; i < X.size(); i++)
    dp += X[i] * Y[i];

  return dp;
} /* End of 'DotProduct' function */

double SqrNorm( std::vector<double> &X ) {
  return DotProduct(X, X);
}

/* Compare two vectors function.
 * ARGUMENTS:
 *   vectors X(N), Y(N);
 * RETURNS: 1 if equal.
 */
int IsEqual( std::vector<double> &X, std::vector<double> &Y, double Epsilon )
{
  int i;
  std::vector<double> a(X.size());

  for (i = 0; i < X.size(); i++)
    a[i] = X[i] - Y[i];

  if (DotProduct(a, a) < Epsilon)
    return 0;
  return (int)(log(DotProduct(a, a)) / log(10.0) + 0.5);
} /* End of 'IsEqual' function */
