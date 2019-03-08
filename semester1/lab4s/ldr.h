#pragma once

#include <vector>

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
               std::vector<double> &x );

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
                       std::vector<std::vector<double>> &R );