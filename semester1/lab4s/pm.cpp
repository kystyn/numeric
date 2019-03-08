#include <math.h>
#include "matr.h"
#include "ldr.h"
#include "pm.h"

/* Pow method algorithm.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> const &A;
 *   - start vector:
 *       std::vector<double> &y0;
 *   - accuracy:
 *       double eps;
 * RETURNS:
 *   (eigen) eigenvalues.
 */
eigen pm( std::vector<std::vector<double>> const &A,
          std::vector<double> const &x0, double eps ) {
  std::vector<double> prev = x0, cur = x0;
  double norm = 0, prevnorm;
  std::vector<std::vector<double>> L, D, R;
  eigen e;

  LDRDecomposition(A, L, D, R);
  printf("ldr: %f\n", matrNorm(MatrAddMatr(MatrMulMatr(L, MatrMulMatr(D, R)), MatrMulNum(A, -1))));

  e.Steps = 0;

  do {
    prev = cur;
    prevnorm = norm;
    LDRSolve(L, D, R, prev, cur);
    norm = VecNormInftySigned(cur);
    cur = VecMulNum(cur, 1.0 / norm);
    e.Steps++;
  } while (fabs(norm - prevnorm) > eps);

  e.Value = norm;
  e.Vector = cur;

  return e;
} /* End of 'pm' function */

/* Pow method with shift algorithm.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> const &A;
 *   - start vector:
 *       std::vector<double> &x0;
 *   - accuracy:
 *       double eps;
 *   - shift:
 *       double shift;
 * RETURNS:
 *   (eigen) eigenvalues.
 */
eigen pmWithShift( std::vector<std::vector<double>> &A, std::vector<double> &x0, double eps, double shift ) {
  std::vector<std::vector<double>> E = A;
  for (int i = 0; i < (int)x0.size(); i++)
    for (int j = 0; j < (int)x0.size(); j++)
      E[i][j] = 0;

  for (int i = 0; i < (int)x0.size(); i++)
    E[i][i] = 1;

  std::vector<std::vector<double>> C = MatrAddMatr(A, MatrMulNum(E, -shift));

  eigen e = pm(C, x0, eps);
  e.Value = e.Value - shift;

  return e;
} /* End of 'pmWithShift' function */