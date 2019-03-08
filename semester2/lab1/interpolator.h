#pragma once

#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include "distribution.h"

class interpolator {
private:
  double (*func)( double );

  std::vector<double (*)( double )> derivatives;
  distribution gridX;
  value_distribution gridY;

  const double epsilon;
  bool
    needRebuildDiff;

  std::vector<std::vector<double>> dividedDiff;

  static unsigned int factorial( unsigned int n ) {
    unsigned int res = 1;
    for (unsigned int i = 1; i <= n; i++)
      res *= i;

    return res;
  }

public:

  interpolator( void ) : func(NULL),
    epsilon(1e-15), needRebuildDiff(true) {}

  interpolator( double (*f)( double ), double a, double b, int nodeCount, 
    std::vector<double (*)( double )> derivatives = std::vector<double (*)( double )>() ) :
    func(f), derivatives(derivatives), epsilon(1e-15), needRebuildDiff(true) {}

  interpolator & operator=( interpolator const &I ) {
    this->derivatives = I.derivatives;
    this->dividedDiff = I.dividedDiff;
    this->func = I.func;
    this->gridX = I.gridX;
    this->gridY = I.gridY;
  }

  interpolator & setFunc( double (*f)( double ), std::vector<double (*)( double )> ders = std::vector<double (*)( double )>() ) {
    func = f;
    derivatives = ders;

    return *this;
  }

  interpolator & setGrid( distribution &newGridX ) {
    gridX = newGridX;
    gridY = value_distribution(gridX, func);
    gridY.eval();

    return *this;
  }

  interpolator & buildDividedDifferenceTable( void ) {
    dividedDiff.resize((derivatives.size() + 1) * gridX.nodeCount);

    for (int i = 0, n = dividedDiff.size(); i < n; i++)
      dividedDiff[i].resize(n - i);

    std::vector<double> extGridX(dividedDiff.size());

    for (int i = 0, n = dividedDiff.size(), k = derivatives.size() + 1; i < n; i++) {
      extGridX[i] = gridX[i / k];
      dividedDiff[0][i] = gridY[i / k];
    }

    for (int i = 1, n = dividedDiff.size(); i < n; i++)
      for (int j = 0, m = dividedDiff[i].size(); j < m; j++) {
        int k = derivatives.size() + 1 - i;

        //if (extGridX[i + j] - extGridX[j] == 0)
        if (k > 0 && j % (derivatives.size() + 1) < k)
          dividedDiff[i][j] = derivatives[i - 1](extGridX[j + i - 1]) / (double)factorial(i);
        else
          dividedDiff[i][j] = (dividedDiff[i - 1][j + 1] - dividedDiff[i - 1][j]) / (extGridX[i + j] - extGridX[j]);
      }
    needRebuildDiff = false;

    return *this;
  }

  double getFunc( double x ) {
    return func(x);
  }

  double HermitPolinom( double x ) {
    if (gridX.needEvalGrid || needRebuildDiff)
      throw "Grid of divided difference table was not rebuilt";

    double res = 0;
    int derCnt = derivatives.size() + 1;

    for (int k = 0, n = dividedDiff.size(); k < n; k++) {
      double localRes = dividedDiff[k][0];

      for (int i = 0; i <= k - 1; i++)
        localRes *= (x - gridX[i / derCnt]);

      res += localRes;
    }

    return res;
  }

  double deviation( void ) {
    if (needRebuildDiff)
      throw "Rebuild divided difference table";

    double d = 0;

    for (int i = 0; i < gridX.nodeCount - 1; i++) {
      double
        //arg = gridX[i],
        arg = (gridX[i] + gridX[i + 1]) / 2.0,
        val = fabs(HermitPolinom(arg) - func(arg));
      if (val > d)
        d = val;
    }

    return d;
  }

  interpolator & output( std::ofstream &f, double step = 1e-2 ) {
    f << (gridX.b - gridX.a) / step << '\n';

    //for (double x = gridX.a; x <= gridX.b; x += step)
    //  f << x << " " << HermitPolinom(x) << "\n";

    f << gridX.nodeCount << " " << deviation() << '\n';
    std::cout << gridX.nodeCount << " " << deviation() << '\n';

    return *this;
  }
};
