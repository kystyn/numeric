#pragma once

#include "def.h"
#include "distribution.h"

class diff_equation_solver {
protected:
  func funct;
  double a, b;
public:
  virtual tabulated_function solve( double tollerance ) {}
};

class euler_cauchy_solver : public diff_equation_solver {
private:
  double cauchyProblem;
  uint frag;

  double deviation( distribution const &valDistr1,
    distribution const &valDistr2 ) const {

    double dev = 0;
    if (valDistr1.getNodeCount() == 2 * valDistr2.getNodeCount()) {
      for (uint i = 0; i < valDistr1.getNodeCount(); i += 2)
        if (fabs(valDistr1[i] - valDistr2[i / 2]) > dev)
          dev = fabs(valDistr1[i] - valDistr2[i / 2]);
    }
    return dev;
  }
public:
  euler_cauchy_solver( void ) {}

  euler_cauchy_solver & setCauchyProblem( double cp ) { cauchyProblem = cp; return *this; }
  euler_cauchy_solver & setFunction( func f ) { funct = f; return *this; }
  euler_cauchy_solver & setBorders( double a, double b ) { this->a = a; this->b = b; return *this; }

  tabulated_function solve( double tollerance ) {
    uint fragmentation = 1;
    uniform argDistr(a, b);
    tabulated_function solution, solutionx2;

    auto eval = [this, fragmentation, &argDistr]( tabulated_function &s ) -> void {
      double h = (b - a) / fragmentation;
      s.setFragmentation(fragmentation);
      s.clear();
      s << cauchyProblem;
      for (uint i = 1; i < fragmentation - 1; i++)
        s << s[i - 1] + h / 2 * (funct(argDistr[i - 1]) + funct(argDistr[i]));

    };

    eval(solutionx2);
    solution << solutionx2[0] + 1;

    for (; deviation(solutionx2, solution) > 3 * tollerance;) {
      solution = solutionx2;
      fragmentation <<= 1;
      eval(solutionx2);
    }

    frag = fragmentation;

    return solutionx2;
  }

  uint getFrag( void ) const { return frag; }
};
