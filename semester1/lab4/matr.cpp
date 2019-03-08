#include "matr.h"

const int matr::OK = 308;

/* Matrix mul vector function.
 * ARGUMENTS:
 *   matrix (H x W);
 *   vector (W);
 *   M, N;
 * RETURNS: vector (M).
 */
vec matr::operator*( vec const &b ) const {
  int i, j;

  if (W != b.getN())
    throw "matr * vec: size mismatch\n";

  vec x(H);

  for (i = 0; i < H; i++) {
    double r = 0;
    for (j = 0; j < W; j++) {
      r += A[i][j] * b[j];
    }
    x[i] = r;
  }

  return x;
}