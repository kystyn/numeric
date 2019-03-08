#include <math.h>
#include "matr.h"
#include "ldr.h"
#include "relax.h"
#include "pm.h"

/* Inverse iterations algorithm.
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
eigen invIt( matr const &A, vec const &x0, double eps, double shift ) {
  eigen e;

  matr L, D, R;
  LDRDecomposition(A, L, D, R);
  printf("ldr: %f\n", !(L * D * R - A));

  e.Steps = 0;
  vec y = x0 / x0.normInftySigned(), x(x0.getN());
  double norm = 0, prevnorm = 1000;

  while (true) {
    x = LDRSolve(L, D, R, y);
    //x = Relax(A, y, x0, 1.2f, eps, n);

    norm = x.normInftySigned();
    y = x / norm;

    if (fabs(1.0 / norm - 1.0 / prevnorm) < eps)
      break;
    printf("||delta_lambda|| = %g\n", fabs(1.0 / norm - 1.0 / prevnorm));

    prevnorm = norm;
    e.Steps++;
  }

  e.Value = 1.0 / norm;
  e.Vector = y;

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
eigen invItWithShift( matr const &A, double eps, double shift ) {
  vec x0(A.getH());

  for (int i = 0; i < x0.getN(); i++)
    x0[i] = 1;

  matr E = A;
  for (int i = 0; i < x0.getN(); i++)
    for (int j = 0; j < x0.getN(); j++)
      E[i][j] = 0;

  for (int i = 0; i < x0.getN(); i++)
    E[i][i] = 1;

  matr C = A - E * shift;

  eigen e = invIt(C, x0, eps, shift);
  e.Value = e.Value + shift;

  return e;
} /* End of 'pmWithShift' function */