#pragma once

#include <vector>
/*
class matr {
private:
  std::vector<std::vector<double>> M;

public:
};*/

/* Matrix mul vector function.
 * ARGUMENTS:
 *   matrix (M x N);
 *   vector (N);
 *   M, N;
 * RETURNS: vector (M).
 */
std::vector<double> MatrMulVec( std::vector<std::vector<double>> const &A, std::vector<double> &b );

/* Matrix mul vector function.
 * ARGUMENTS:
 *   matrix (M x N);
 *   vector (N);
 *   M, N;
 * RETURNS: vector (M).
 */
std::vector<double> VecMulMatr( std::vector<double> &b, std::vector<std::vector<double>> const &A );

/* Vec mul num function.
 */
std::vector<double> VecMulNum( std::vector<double> &x, double Num );

/* Vec add vec  function.
 */
std::vector<double> VecAddVec( std::vector<double> &x, std::vector<double> &y );

/* Matrix mul matrix function.
 * ARGUMENTS:
 *   matrix (M x N), (N X K);
 *   M, N, K;
 * RETURNS: matrix (M X K).
 */
std::vector<std::vector<double>> MatrMulMatr( std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> const &B );

/* Matrix mul matrix function.
 * ARGUMENTS:
 *   matrix (M x N), (N X K);
 *   M, N, K;
 * RETURNS: matrix (M X K).
 */
std::vector<std::vector<double>> MatrMulNum( std::vector<std::vector<double>> const &A, double M );

/* Transpose matrix function.
 * ARGUMENTS:
 *   - matrix M X N;  
 *   - M, N;
 * RETURNS: matrix N X M.
 */
std::vector<std::vector<double>> Transposing( std::vector<std::vector<double>> const &A );

/* Dot product function.
 * ARGUMENTS:
 *   - vectors:
 *       X, Y (N);
 * RETURNS: dot product.
 */
double DotProduct( std::vector<double> &X, std::vector<double> &Y );

double SqrNorm( std::vector<double> &X );

double VecNormInftySigned( std::vector<double> &X );

double matrNorm( std::vector<std::vector<double>> const &A );

/* Compare two vectors function.
 * ARGUMENTS:
 *   vectors X(N), Y(N);
 * RETURNS: 1 if equal.
 */
int IsEqual( std::vector<double> &X, std::vector<double> &Y, double Epsilon );

std::vector<std::vector<double>> MatrAddMatr( std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B );
