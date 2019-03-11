#pragma once

#include <vector>
#include <cmath>

#include "def.h"

/* Matrix class */
namespace mth {
class vec;
static const double Tollerance = 1e-15;
class matr {
private:
  std::vector<std::vector<double>> A;
  size_t H, W;

public:
  static const int OK;
  matr() : H(0), W(0) {}

  matr( std::vector<std::vector<double>> const &A ) : A(A), H(A.size()), W(A[0].size()) {}

  matr( matr const &B ) {
    A = B.A;
    H = B.H;
    W = B.W;
  }

  matr & operator=( matr const &B ) {
    A = B.A;
    H = B.H;
    W = B.W;

    return *this;
  }

  matr( size_t H, size_t W = 0 ) : H(H), W(W) {
    A.resize(H);

    if (W == 0)
        this->W = H;

    for (unsigned int i = 0; i < H; i++)
      A[i].resize(this->W);
  }

  size_t getW() const { return W; }

  size_t getH() const { return H; }


  std::vector<double> & operator[]( unsigned int i ) {
    return A[i];
  }

  std::vector<double> operator[]( unsigned int i ) const {
    return A[i];
  }

  /* Matrix mul vector function.
   * ARGUMENTS:
   *   matrix (H x W);
   *   vector (W);
   *   M, N;
   * RETURNS: vector (M).
   */
  vec operator*( vec const &b ) const;

  /* Matrix mul matrix function.
   * ARGUMENTS:
   *   matrix (M x N), (N X K);
   *   M, N, K;
   * RETURNS: matrix (M X K).
   */
  matr operator*( matr const &B ) const {
    if (W != B.H)
      throw "matr * matr: size mismatch\n";

    matr C(H, B.W);

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < B.W; j++)
      {
        C[i][j] = 0;
        for (unsigned int k = 0; k < W; k++)
          C[i][j] += A[i][k] * B[k][j];
      }

    return C;
  }

  matr operator+( matr const &B ) const {
    if (H != B.H && W != B.W)
      throw "matr + matr: size mismatch\n";

    matr C(H, W);
    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        C[i][j] = A[i][j] + B[i][j];

    return C;
  }

  matr & operator+=( matr const &B ) {
    if (H != B.H && W != B.W)
      throw "matr + matr: size mismatch\n";

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        A[i][j] += B[i][j];

    return *this;
  }

  matr operator-( matr const &B ) const {
    if (H != B.H && W != B.W)
      throw "matr + matr: size mismatch\n";

    matr C(H, W);
    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        C[i][j] = A[i][j] - B[i][j];

    return C;
  }

  matr & operator-=( matr const &B ) {
    if (H != B.H && W != B.W)
      throw "matr + matr: size mismatch\n";

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        A[i][j] -= B[i][j];

    return *this;
  }

  /* Matrix mul matrix function.
   * ARGUMENTS:
   *   matrix (M x N), (N X K);
   *   M, N, K;
   * RETURNS: matrix (M X K).
   */
  matr operator*( double Num ) const {
    matr C = A;

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        C[i][j] *= Num;

    return C;
  }

  /* Matrix mul matrix function.
   * ARGUMENTS:
   *   matrix (M x N), (N X K);
   *   M, N, K;
   * RETURNS: matrix (M X K).
   */
  matr & operator*=( double Num ) {

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        A[i][j] *= Num;

    return *this;
  }

  matr transposing( void ) const {
    matr At(W, H);

    for (unsigned int i = 0; i < W; i++)
      for (unsigned int j = 0; j < H; j++)
        At[i][j] = A[j][i];

    return At;
  }

  matr & transpose( void ) {
    matr At(W, H);

    for (unsigned int i = 0; i < W; i++)
      for (unsigned int j = 0; j < H; j++)
        At[i][j] = A[j][i];

    *this = At;

    return *this;
  }

  double operator!( void ) const {
    double n = 0;

    for (unsigned int i = 0; i < H; i++)
      for (unsigned int j = 0; j < W; j++)
        n += A[i][j] * A[i][j];
    
    return sqrt(n);
  }
};

class vec {
private:
  std::vector<double> V;
  size_t N;
public:
  vec() : N(0) {}

  vec( std::vector<double> const &V ) : V(V), N(V.size()) {}

  vec( size_t N ) : V(N), N(N) {}

  vec( size_t N, int Val ) : V(N, Val), N(N) {}

  size_t getN() const { return N; }

  double & operator[]( unsigned int i ) {
    return V[i];
  }

  double operator[]( unsigned int i ) const {
    return V[i];
  }

  /* Vector mul matrix function.
   * ARGUMENTS:
   *   vector (transposed) (M);
   *   matrix (M x N);
   *   M, N;
   * RETURNS: vector (transposed) (N).
   */
  vec operator*( matr const &A ) const {
    if (N != A.getH())
      throw "vec * matr: size mismatch\n";

    vec x(A.getW());
    
    for (unsigned int i = 0; i < N; i++)
    {
      double r = 0;
      for (unsigned int j = 0; j < A.getW(); j++)
        r += V[i] * A[j][i];
      x[i] = r;
    }

    return x;
  }

  /* Vector mul num function.
   * ARGUMENTS:
   *   vector (transposed) (M);
   *   number;
   * RETURNS: vector.
   */
  vec operator*( double Num ) const {
    vec x = V;
    
    for (unsigned int i = 0; i < N; i++)
      x[i] *= Num;

    return x;
  }

  /* Vector mul num function.
   * ARGUMENTS:
   *   vector (transposed) (M);
   *   number;
   * RETURNS: vector.
   */
  vec & operator*=( double Num ) {
    for (unsigned int i = 0; i < N; i++)
      V[i] *= Num;

    return *this;
  }

  /* Vector mul num function.
   * ARGUMENTS:
   *   vector (transposed) (M);
   *   number;
   * RETURNS: vector.
   */
  vec operator/( double Num ) const {
    return *this * (1.0 / Num);
  }

  /* Vector mul num function.
   * ARGUMENTS:
   *   vector (transposed) (M);
   *   number;
   * RETURNS: vector.
   */
  vec & operator/=( double Num ) {
    *this *= 1.0 / Num;

    return *this;
  }

  vec operator+( vec const &b ) const {
    if (b.getN() != N)
      throw "vec + vec: size mismatch!\n";

    vec a = V;

    for (unsigned int i = 0; i < N; i++)
      a[i] += b[i];

    return a;
  }

  vec operator+=( vec const &b ) {
    if (b.getN() != N)
      throw "vec + vec: size mismatch!\n";

    for (unsigned int i = 0; i < N; i++)
      V[i] += b[i];

    return *this;
  }

  vec operator-( vec const &b ) const {
    if (b.getN() != N)
      throw "vec + vec: size mismatch!\n";

    vec a = V;

    for (unsigned int i = 0; i < N; i++)
      a[i] -= b[i];

    return a;
  }

  vec operator-=( vec const &b ) {
    if (b.getN() != N)
      throw "vec + vec: size mismatch!\n";

    for (unsigned int i = 0; i < N; i++)
      V[i] -= b[i];

    return *this;
  }

  double operator*( vec const &x ) const {
    if (x.getN() != N)
      throw "vec * vec: size mismatch!\n";

    double r = 0;

    for (unsigned int i = 0; i < N; i++)
      r += V[i] * x[i];

    return r;
  }

  double operator!( void ) const {
    return sqrt(*this * *this);
  }

  double normInftySigned( void ) const {
    double t = 0;

    for (unsigned int i = 0; i < N; i++)
      if (fabs(V[i]) > fabs(t))
        t = V[i];

    return t;
  }
  double normInfty( void ) const {
    return fabs(normInftySigned());
  }

  /* Compare two vectors function.
   * ARGUMENTS:
   *   vectors X(N), Y(N);
   * RETURNS: 0 if equal.
   */
  int cmp( vec const &x, double epsilon ) const {
    vec a = *this - x;
    double sqrnorm = a * a;
    if (sqrnorm < epsilon)
      return matr::OK;
    return (int)(log(sqrnorm) / log(10.0) + 0.5);
  }

  operator std::vector<double>( void ) {
    return V;
  }
}; /* End of 'vec' class */
} // end of 'mth' namespace
