#include "ldr.h"
#include "simple.h"
#include "matr.h"

/* LDR decomposition function.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> &A;
 *   - size:
 *       int N;
 *   - right vector:
 *       std::vector<double> &b;
 *   - non-allocated L, D, R matrices;
 * RETURNS: None.
 */
void LDRDecomposition( std::vector<std::vector<double>> const &A,
                       std::vector<std::vector<double>> &L,
                       std::vector<std::vector<double>> &D,
                       std::vector<std::vector<double>> &R )
{
  int i, k, N = A.size();

  /* Variables for LDR */
  std::vector<double> f(N), g(N), x(N), y(N); 
  double beta;

  L.resize(N);
  D.resize(N);
  R.resize(N);

  for (i = 0; i < N; i++) {
    L[i] = std::vector<double>(N);
    D[i] = std::vector<double>(N);
    R[i] = std::vector<double>(N);
  }
  
  /* A1 = L1 * D1 * R1 */

  L[0][0] = 1;
  D[0][0] = A[0][0];
  R[0][0] = 1;

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
    Left(MatrMulMatr(Transposing(R), Transposing(D)),
      g, x, k - 1);

    for (i = 0; i < k - 1; i++)
    {
      L[k - 1][i] = x[i];
      L[i][k - 1] = 0;
    }
    L[k - 1][k - 1] = 1;

    /*** Set Rk ***/
    /* Solve Y */
    Left(MatrMulMatr(L, D), f, y, k - 1);

    for (i = 0; i < k - 1; i++)
    {
      R[i][k - 1] = y[i];
      R[k - 1][i] = 0;
    }
    R[k - 1][k - 1] = 1;

    /*** Set Dk ***/
    /* Solve beta */
    beta = A[k - 1][k - 1] - DotProduct(VecMulMatr(x, D), y);
    for (i = 0; i < k - 1; i++)
    {
      D[i][k - 1] = 0;
      D[k - 1][i] = 0;
    }
    D[k - 1][k - 1] = beta;
  }
} /* End of 'LDRDecomposition' function */

/* LDR solve function.
 * ARGUMENTS:
 *   - matrix (N X N);
 *   - right vector b (N);
 *   - N;
 *   - allocated result vector x (N);
 * RETURNS: None.
 */
void LDRSolve( std::vector<std::vector<double>> &L,
               std::vector<std::vector<double>> &D,
               std::vector<std::vector<double>> &R,
               std::vector<double> &b,
               std::vector<double> &x )
{
  int N = L.size();
  std::vector<double> y(N), z(N);

  /* Ly = b */
  Left(L, b, y, N);
  /* Dz = y */
  Diagonal(D, y, z, N);
  /* Rx = z */
  Right(R, z, x, N);
} /* End of 'LDRSolve' function */