#pragma once

#include <vector>

namespace mth {
class matr;
class vec;

namespace lieqsys {

/* LDR solve function.
 * ARGUMENTS:
 *   - matrix (N X N);
 *   - right vector b (N);
 *   - N;
 *   - allocated result vector x (N);
 * RETURNS: None.
 */
vec LDRSolve( matr const &L,
               matr const &D,
               matr const &R,
               vec const &b );

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
void LDLTDecomposition( matr const &A,
                       matr &L,
                       matr &D,
                       matr &R );
}
}
