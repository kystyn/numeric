#pragma once

#include <map>
#include <fstream>
#include <string>
#include "distribution.h"
#include "interpolator.h"

struct data {
  double a, b;
  int nodeCount;
  char gridT;
  int derivativeT;
  int funcT;

  data( void ) : a(0), b(0), nodeCount(0), gridT('u'), derivativeT(0), funcT(0){}

  data( double a, double b, int node, char gr, int der, int f ) :
  a(a), b(b), nodeCount(node), gridT(gr), derivativeT(der), funcT(f) {}
};

class controller {
private:
  interpolator i;

  std::map<char, distribution *> availableGridX;
  std::vector<double (*)( double )> availableFunctions;
  std::vector<std::vector<double (*)( double )>> availableDerivatives;
  std::vector<data> loadedData;

  controller & loadFromFile( const char *fileName ) {
    std::ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      double a, b;
      char gridT;
      int useder, funcT;
      int nc;
      if (f >> a)
        f >> b >> nc >> gridT >> useder >> funcT;
      else
        break;
      loadedData.push_back(data(a, b, nc, gridT, useder, funcT));
    }

    return *this;
  }

public:
  static double diffFunc( double x ) {
    return log(1 + cos(x) * cos(x));
    //return pow(x, 8) + 1;
  }

  static double nonDiffFunc( double x ) {
    return fabs(log(1 + cos(x) * cos(x)) - 0.1);
    //return pow(x, 8) + 1;
  }

  static double derivativeOfDiffFunc( double x ) {
    return -sin(2 * x) / (1 + cos(x) * cos(x));
    //return 8 * pow(x, 7);
  }

  static double pow_1( double a ) {
    return a == 0 ? 1 : -1;
  }

  static double derivativeOfNonDiffFunc( double x ) {
    return (double)pow_1(x < 3) * -sin(2 * x) / (1 + cos(x) * cos(x));
    //return 8 * pow(x, 7);
  }

  distribution *& operator[]( char p ) {
    return availableGridX[p];
  }

  controller & operator<<( double (*v)( double ) ) {
    availableFunctions.push_back(v);
    return *this;
  }

  controller & operator<<( std::vector<double (*)( double )> &v ) {
    availableDerivatives.push_back(v);
    return *this;
  }

  controller( const char *fileName ) {
    loadFromFile(fileName);
  }

  void run( const char *fileName ) {
    auto fs = std::ofstream(fileName);

    for (auto &d : loadedData) {
      auto f = availableFunctions[d.funcT];
      auto df = availableDerivatives[d.derivativeT];
      i.setFunc(f, df);

      i.setGrid(availableGridX[d.gridT]->setBorders(d.a, d.b, d.nodeCount).eval());
      i.buildDividedDifferenceTable();
      fs << d.gridT << ' ' << d.funcT << ' ' << d.derivativeT << '\n';
      i.output(fs);
    }
  }
};
