#pragma once

#include "def.h"
#include "distribution.h"

class diff_equation_solver {
protected:
  func2var funct;
  double a, b;
  uniform interval;
public:
  virtual tabulated_function solve( double tollerance ) = 0;
};

class euler_cauchy_solver : public diff_equation_solver {
private:
  double cauchyProblem;
  uint minFrag;
  uint maxFrag;

  double deviation( distribution const &valDistr1,
    distribution const &valDistr2 ) const {

    return fabs(valDistr1[valDistr1.getNodeCount() - 1] - valDistr2[valDistr2.getNodeCount() - 1]);
  }
public:
  euler_cauchy_solver( void ) {}

  euler_cauchy_solver & setCauchyProblem( double cp ) { cauchyProblem = cp; return *this; }
  euler_cauchy_solver & setFragmentation( uint frag ) {
      interval.setBorders(frag);
      interval.eval();
      return *this;
  }

  euler_cauchy_solver & setFunction( func2var f ) { funct = f; return *this; }
  euler_cauchy_solver & setBorders( double a, double b ) {
      this->a = a;
      this->b = b;
      interval.setBorders(a, b);
      interval.eval();
      return *this;
  }

  tabulated_function solve( double tollerance ) {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;

    solution << make_pair<>(a, cauchyProblem);

    auto eval = [this](
            tabulated_function &s, double x1, double cp, double x2, uint localFrag ) -> void {
      uniform argDistr(x1, x2);
      double h = (x2 - x1) / (localFrag - 1);
      argDistr.setBorders(localFrag);
      argDistr.eval();
      s.clear();

      s << make_pair<>(x1, cp);

      for (uint i = 1; i < localFrag; i++)
        s << make_pair<>(
                 argDistr[i],
                 s[i - 1] + h / 2 * (funct(argDistr[i - 1], s[i - 1]) + funct(argDistr[i], s[i - 1] + h * funct(argDistr[i - 1], s[i - 1]))));

    };

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + 1);
        eval(s2, interval[i], solution[i], interval[i + 1], localFrag);

        for (; deviation(s2, s1) > 3 * tollerance;) {
            s1 = s2;
            eval(s2, interval[i], solution[i], interval[i + 1], localFrag <<= 1);
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }

  uint getMinFrag( void ) const { return minFrag; }
  uint getMaxFrag( void ) const { return maxFrag; }
};
