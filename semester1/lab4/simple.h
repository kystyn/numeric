#pragma once
#include "matr.h"

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
void Diagonal( matr const &A, vec const &b, vec &x, int N );

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
void Left( matr const &A, vec const &b, vec &x, int N );

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
void Right( matr const &A, vec const &b, vec &x, int N );
