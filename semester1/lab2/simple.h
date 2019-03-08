#pragma once

/* Diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Diagonal( double *A[], double *b, double *x, int N );

/* Left diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Left( double **A, double *b, double *x, int N );

/* Right diagonal system solve function.
 * ARGUMENTS:
 *   - matrix:
 *       double **A;
 *   - size:
 *       int N;
 *   - right vector:
 *       double *b;
 *   - result (allocated) vector:
 *       double *x;
 * RETURNS: None.
 */
void Right( double *A[], double *b, double *x, int N );
