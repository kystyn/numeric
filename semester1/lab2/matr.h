#pragma once

/* Matrix mul vector function.
 * ARGUMENTS:
 *   matrix (M x N);
 *   vector (N);
 *   M, N;
 * RETURNS: vector (M).
 */
double * MatrMulVec( double **A, double *b, int M, int N );

/* Vector mul matrix function.
 * ARGUMENTS:
 *   vector (transposed) (M);
 *   matrix (M x N);
 *   M, N;
 * RETURNS: vector (transposed) (N).
 */
double * VecMulMatr( double *b, double **A, int M, int N );

/* Vec mul num function.
 */
double * VecMulNum( double *x, int N, double Num );

/* Matrix mul matrix function.
 * ARGUMENTS:
 *   matrix (M x N), (N X K);
 *   M, N, K;
 * RETURNS: matrix (M X K).
 */
double ** MatrMulMatr( double **A, double **B, int M, int N, int K );

/* Transpose matrix function.
 * ARGUMENTS:
 *   - matrix M X N;  
 *   - M, N;
 * RETURNS: matrix N X M.
 */
double ** Transposing( double **A, int M, int N );

/* Dot product function.
 * ARGUMENTS:
 *   - vectors:
 *       X, Y (N);
 * RETURNS: dot product.
 */
double DotProduct( double *X, double *Y, int N );

/* Compare two vectors function.
 * ARGUMENTS:
 *   vectors X(N), Y(N);
 * RETURNS: 1 if equal.
 */
int IsEqual( double *X, double *Y, int N );
