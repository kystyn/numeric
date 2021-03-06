#pragma once

#include <fstream>
#include <string>
#include <random>
#include <memory>
#include <vector>

#include "distribution.h"
#include "integrator.h"

namespace kystyn {
struct data {
  double a, b;
  uint funcT;
  double Tollerance;

  data( void ) : a(0), b(0), funcT(0), Tollerance(1e-16) {}

  data( double a, double b, uint funcT, double Tollerance ) :
    a(a), b(b),
    funcT(funcT), Tollerance(Tollerance) {}
};
}

class controller {
private:
  rado_integral_n2 i;

  vector<func> availableFunctions;
  vector<kystyn::data> loadedData;

  controller & loadFromFile( const char *fileName ) {
    ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      double a, b;
      uint funcT;
      double tollerance;

      if (f >> a)
        f >> b >> funcT >> tollerance;
      else
        break;
      loadedData.push_back(kystyn::data(a, b, funcT, pow(10, tollerance)));
    }

    return *this;
  }

public:

  static double diffFunc( double x ) {
    return x * x * x * cos(x);
    //return pow(x, 8) + 1;
  }

  static double nonDiffFunc( double x ) {
    return fabs(x * x * x * cos(x));
  }

  controller( const char *fileName = nullptr ) {
      if (fileName)
        loadFromFile(fileName);
  }

  controller & operator<<( func v ) {
    availableFunctions.push_back(v);
    return *this;
  }

  void run( const char *fileName ) {
    auto fs = ofstream(fileName);

    for (auto &d : loadedData) {
        auto f = availableFunctions[d.funcT];
        fs << std::scientific << i.
            setBorders(d.a, d.b).
            setTollerance(d.Tollerance).setFunction(f)() << ' ' << i.getFragmentation() << std::endl;
    }
  }
};
