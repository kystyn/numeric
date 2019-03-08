#pragma once
#include <vector>

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
void Diagonal( std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x, int N );

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
void Left( std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x, int N );

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
void Right( std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x, int N );
