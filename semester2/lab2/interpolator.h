#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>

#include "distribution.h"
#include "tabulated_basis_func.h"

#include "matr.h"
#include "ldr.h"
#include "relax.h"

class least_square_interpolator {
public:
    using basis_func = std::function<double(double, int)>;
private:
  func Func;
  basis_func BasisFunc;

  distribution GridX;
  value_distribution GridY;

  distribution Weights;

  std::vector<double> PolinomialCoeffs;
  size_t PolinomDimension;

  static unsigned int factorial( unsigned int n ) {
    unsigned int res = 1;
    for (unsigned int i = 1; i <= n; i++)
      res *= i;

    return res;
  }

public:

  least_square_interpolator( void ) : Func(nullptr) {}

  least_square_interpolator( func F, basis_func BF, size_t PolinomDimension, distribution const &GridX ) :
      Func(F), BasisFunc(BF), GridX(GridX), PolinomDimension(PolinomDimension + 1) {
      GridY = value_distribution(GridX, Func);
  }

  least_square_interpolator & operator=( least_square_interpolator const &I ) {
    this->Func = I.Func;
    this->GridX = I.GridX;
    this->GridY = I.GridY;

    return *this;
  }

  least_square_interpolator & setFunc( func F ) {
    Func = F;

    return *this;
  }

  least_square_interpolator & setBasisFunc( basis_func F, size_t PolinomDim = 0 ) {
    BasisFunc = F;
    PolinomDimension = PolinomDim;

    return *this;
  }
  least_square_interpolator & setPolyDimension( size_t PolinomDim ) {
      PolinomDimension = PolinomDim;

      return *this;
  }

  least_square_interpolator & setGrid( distribution &NewGridX ) {
    GridX = NewGridX;
    GridY = value_distribution(GridX, Func);

    return *this;
  }

  least_square_interpolator & setWeights( std::vector<double> const &NWeights ) {
      Weights = NWeights;

      return *this;
  }

  void genPolinom( void ) {
      if (PolinomDimension > GridX.NodeCount)
          throw "There should be less weights than grid nodes";

      if (Weights.NodeCount != GridX.NodeCount)
          throw "There should be the same quantity of weights as grid size";

      mth::matr A(PolinomDimension);
      mth::vec b(PolinomDimension);

      tabulated_basis_func y(GridY, Weights);

      for (uint i = 0; i < A.getH(); i++) {
          tabulated_basis_func
            bi(BasisFunc, GridX, i, Weights);
          for (uint j = 0; j < A.getW(); j++) {
              tabulated_basis_func bj(BasisFunc, GridX, j, Weights);
              A[i][j] = bi * bj;
          }
          b[i] = y * bi;
      }
      mth::matr L, D, R;
      mth::lieqsys::LDLTDecomposition(A, L, D, R);
      cout << "norm " << !(A - L * D * R) << endl;
      PolinomialCoeffs = mth::lieqsys::LDRSolve(L, D, R, b);
      ///mth::lieqsys::Relax(A, b, mth::vec(b.getN(), 0), 1.2, Tolerance, steps);
  }

  double operator()( double X ) {
      if (PolinomialCoeffs.size() == 0)
          throw "Can not evaluate polinom!";

      double res = 0;
      double PoweredX = 1;

      for (auto const &c : PolinomialCoeffs) {
          res += PoweredX * c;
          PoweredX *= X;
      }

      return res;
  }

  double operator!( void ) {
      double m = 0;
      for (uint i = 0; i < GridX.NodeCount - 1; i++) {
          double arg = (GridX[i] + GridX[i + 1]) / 2.0;
          m = std::max(m, fabs((*this)(arg) - Func(arg)));
      }

      return m;
  }

  least_square_interpolator & output( std::ofstream &f ) {

    f << GridX.NodeCount << " " << !*this << '\n';
    //std::cout << GridX.NodeCount << " " << !*this << '\n';

    return *this;
  }
};
