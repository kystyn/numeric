#pragma once

#include <vector>

struct eigen {
  double
    Value;
  std::vector<double> Vector;
  int Steps;
};

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
eigen pmWithShift( std::vector<std::vector<double>> &A, std::vector<double> &x0, double eps, double shift );