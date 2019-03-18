#include "ldr.h"
#include "simple.h"
#include "matr.h"

using namespace mth;

/* LDR decomposition function.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> &A;
 *   - size:
 *       uint N;
 *   - right vector:
 *       std::vector<double> &b;
 *   - non-allocated L, D, R matrices;
 * RETURNS: None.
 */
void lieqsys::LDLTDecomposition( matr const &A,
                       matr &L,
                       matr &D,
                       matr &R )
{
  size_t N = A.getH();

  /* Variables for LDR */
  vec g(N), x(N);
  double beta;

  L = matr(N, N);
  D = matr(N, N);

  /* A1 = L1 * D1 * R1 */

  L[0][0] = 1;
  D[0][0] = A[0][0];

  for (uint k = 2; k <= N; k++)
  {
    /* Fill g_k-1T */
    for (uint i = 0; i < k - 1; i++)
      g[i] = A[k - 1][i];

    /*** Set Lk ***/
    /* Solve X */
    mth::lieqsys::Left(L * D,
      g, x, k - 1);

    for (uint i = 0; i < k - 1; i++)
    {
      L[k - 1][i] = x[i];
      L[i][k - 1] = 0;
    }
    L[k - 1][k - 1] = 1;

    /*** Set Dk ***/
    /* Solve beta */
    beta = A[k - 1][k - 1] - x * D * x;
    for (uint i = 0; i < k - 1; i++)
    {
      D[i][k - 1] = 0;
      D[k - 1][i] = 0;
    }
    D[k - 1][k - 1] = beta;
  }
  R = L.transposing();
} /* End of 'LDRDecomposition' function */

/* LDR solve function.
 * ARGUMENTS:
 *   - matrix (N X N);
 *   - right vector b (N);
 *   - N;
 *   - allocated result vector x (N);
 * RETURNS: None.
 */
vec lieqsys::LDRSolve( matr const &L,
               matr const &D,
               matr const &R,
               vec const &b )
{
  size_t N = L.getH();
  vec x(N), y(N), z(N);

  /* Ly = b */
  mth::lieqsys::Left(L, b, y, b.getN());
  /* Dz = y */
  mth::lieqsys::Diagonal(D, y, z, b.getN());
  /* Rx = z */
  mth::lieqsys::Right(R, z, x, b.getN());

  return x;
} /* End of 'LDRSolve' function */
