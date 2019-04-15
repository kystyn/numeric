#pragma once

#include "def.h"
#include "distribution.h"

class diff_equation_solver {
protected:
  func2var funct;
  double a, b;
public:
  virtual tabulated_function solve( double tollerance ) = 0;
};

class euler_cauchy_solver : public diff_equation_solver {
private:
  double cauchyProblem;
  uint minFrag;
  uint maxFrag;

  double maxDeviation( distribution const &valDistr1,
    distribution const &valDistr2 ) const {

    double dev = 0;
    for (uint i = 0; i < valDistr1.getNodeCount(); i += 2)
      if (fabs((long double)(valDistr1[i] - valDistr2[i / 2])) > dev)
        dev = fabs((long double)(valDistr1[i] - valDistr2[i / 2]));

    return dev;
  }

  bool reachedMinDeviation( distribution const &valDistr1,
    distribution const &valDistr2, double tollerance ) const {

    for (uint i = 2; i < valDistr1.getNodeCount(); i += 2)
      if (fabs((long double)(valDistr1[i] - valDistr2[i / 2])) <= 3 * tollerance)
        return true;

    return false;
  }
public:
  euler_cauchy_solver( void ) {}

  euler_cauchy_solver & setCauchyProblem( double cp ) { cauchyProblem = cp; return *this; }
  euler_cauchy_solver & setFunction( func2var f ) { funct = f; return *this; }
  euler_cauchy_solver & setBorders( double a, double b ) { this->a = a; this->b = b; return *this; }

  tabulated_function solve( double tollerance ) {
    uint fragmentation = 2;
    minFrag = fragmentation;
    uniform argDistr(a, b);
    tabulated_function solution, solutionx2;

    solution.setArgumentDistrubution(shared_ptr<distribution>(new uniform(a, b)));
    solutionx2.setArgumentDistrubution(shared_ptr<distribution>(new uniform(a, b)));

    auto eval = [this, &fragmentation, &argDistr]( tabulated_function &s ) -> void {
      double h = (b - a) / (fragmentation - 1);
      s.setFragmentation(fragmentation);
      argDistr.setBorders(fragmentation);
      argDistr.eval();
      s.clear();

      s << cauchyProblem;

      for (uint i = 1; i < fragmentation; i++)
        s << s[i - 1] + h / 2 * (funct(argDistr[i - 1], s[i - 1]) + funct(argDistr[i], s[i - 1] + h * funct(argDistr[i - 1], s[i - 1])));

    };

    eval(solutionx2);
    solution << solutionx2[0] + 1;

    bool evalMinDev = true;

    for (; maxDeviation(solutionx2, solution) > 33 * tollerance;) {
      if (evalMinDev)
        if (reachedMinDeviation(solutionx2, solution, tollerance)) {
          evalMinDev = false;
          minFrag = fragmentation;
        }
      solution = solutionx2;
      fragmentation <<= 1;
      eval(solutionx2);
    }

    maxFrag = fragmentation;

    return solutionx2;
  }

  uint getMinFrag( void ) const { return minFrag; }
  uint getMaxFrag( void ) const { return maxFrag; }
};
