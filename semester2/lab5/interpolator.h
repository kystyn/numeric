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
  func f;

  std::vector<double (*)( double )> derivatives;
  distribution gridX;
  grid gridY;

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

  interpolator( void ) : f(nullptr),
    needRebuildDiff(true) {}

  interpolator( double (*f)( double ),
    std::vector<double (*)( double )> derivatives = std::vector<double (*)( double )>() ) :
    f(f), derivatives(derivatives), needRebuildDiff(true) {}

  interpolator & operator=( interpolator const &I ) {
    this->derivatives = I.derivatives;
    this->dividedDiff = I.dividedDiff;
    this->f = I.f;
    this->gridY = I.gridY;

    return *this;
  }

  interpolator & setFunc( double (*f)( double ), std::vector<double (*)( double )> ders = std::vector<double (*)( double )>() ) {
    f = f;
    derivatives = ders;

    return *this;
  }

  interpolator & setGridX( distribution const &newGrid ) {
    gridX = newGrid;

    return *this;
  }

  interpolator & setGridY( grid const &newGrid ) {
    gridY = newGrid;

    return *this;
  }

  interpolator & buildDividedDifferenceTable( void ) {
    dividedDiff.resize((derivatives.size() + 1) * gridY.getNodeCount());

    for (long long i = 0, n = dividedDiff.size(); i < n; i++)
      dividedDiff[i].resize(n - i);

    std::vector<double> extGridX(dividedDiff.size());

    for (long long i = 0, n = dividedDiff.size(), k = derivatives.size() + 1; i < n; i++) {
      extGridX[i] = gridX[static_cast<unsigned int>(i / k)];
      dividedDiff[0][i] = gridY[static_cast<unsigned int>(i / k)];
    }

    for (long long i = 1, n = dividedDiff.size(); i < n; i++)
      for (long long j = 0, m = dividedDiff[i].size(); j < m; j++) {
        long long k = static_cast<long long>(derivatives.size() + 1 - i);

        //if (extGridX[i + j] - extGridX[j] == 0)
        if (k > 0 && j % (derivatives.size() + 1) < k)
          dividedDiff[i][j] = derivatives[i - 1](extGridX[j + i - 1]) / double(factorial(i));
        else
          dividedDiff[i][j] = (dividedDiff[i - 1][j + 1] - dividedDiff[i - 1][j]) / (extGridX[i + j] - extGridX[j]);
      }
    needRebuildDiff = false;

    return *this;
  }

  double getFunc( double x ) {
    return f(x);
  }

  func HermitPolinom( void ) {
    auto f = [this]( double x ) {

      double res = 0;
      long long derCnt = derivatives.size() + 1;

      for (long long k = 0, n = dividedDiff.size(); k < n; k++) {
        double localRes = dividedDiff[k][0];

        for (long long i = 0; i <= k - 1; i++)
          localRes *= (x - gridX[i / derCnt]);

        res += localRes;
      }

      return res;
    };
    return f;
  }

  double deviation( void ) {
    if (needRebuildDiff)
      throw "Rebuild divided difference table";

    double d = 0;

    for (long long i = 0; i < gridX.getNodeCount() - 1; i++) {
      double
        //arg = gridX[i],
        arg = (gridX[i] + gridX[i + 1]) / 2.0,
        val = fabs(HermitPolinom()(arg) - f(arg));
      if (val > d)
        d = val;
    }

    return d;
  }

  interpolator & output( std::ofstream &f, double step = 1e-2 ) {
    f << (gridX.getB() - gridX.getA()) / step << '\n';

    //for (double x = gridX.a; x <= gridX.b; x += step)
    //  f << x << " " << HermitPolinom(x) << "\n";

    f << gridX.getNodeCount() << " " << deviation() << '\n';
    std::cout << gridX.getNodeCount() << " " << deviation() << '\n';

    return *this;
  }
};
