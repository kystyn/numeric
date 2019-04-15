#pragma once

#include <fstream>
#include <string>
#include <random>
#include <memory>
#include <vector>
#include <iomanip>

#include "distribution.h"
#include "diff_equations.h"

namespace kyst {
struct data {
  double a, b;
  double cauchyProblem;
  double Tollerance;

  data( void ) : a(0), b(0), cauchyProblem(0), Tollerance(1e-16) {}

  data( double a, double b, double cauchyProblem, double Tollerance ) :
    a(a), b(b),
    cauchyProblem(cauchyProblem), Tollerance(Tollerance) {}
};
}

class controller {
private:
  euler_cauchy_solver ecs;
  func2var funct;
  vector<kyst::data> loadedData;

  controller & loadFromFile( const char *fileName ) {
    ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      double a, b;
      double cp;
      double tollerance;

      if (f >> a)
        f >> b >> cp >> tollerance;
      else
        break;
      loadedData.push_back(kyst::data(a, b, cp, pow(10, tollerance)));
    }

    return *this;
  }

public:

  static double diffFunc( double x ) {
    return sin(x);
    //return pow(x, 8) + 1;
  }

  static double nonDiffFunc( double x ) {
    return fabs(sin(x));
  }

  controller( const char *fileName = nullptr ) {
      if (fileName)
        loadFromFile(fileName);
  }

  controller & operator<<( func2var v ) {
    funct = v;
    return *this;
  }

  void run( const char *fileName ) {
    auto fs = ofstream(fileName);

    for (auto &d : loadedData) {
        fs << std::setprecision(16) << ecs.
            setBorders(d.a, d.b).
            setFunction(funct).setCauchyProblem(d.cauchyProblem).
            solve(d.Tollerance) << ecs.getFrag() << std::endl;
    }
  }
};
