#pragma once

#include <vector>
#include "matr.h"

struct eigen {
  double
    Value;
  vec Vector;
  int Steps;
};

/* Pow method with shift algorithm.
 * ARGUMENTS:
 *   - matrix:
 *       std::vector<std::vector<double>> const &A;
 *   - accuracy:
 *       double eps;
 *   - shift:
 *       double shift;
 * RETURNS:
 *   (eigen) eigenvalues.
 */
eigen invItWithShift( matr const &A, double eps, double shift );