#pragma once

#include "def.h"
#include "distribution.h"

class diff_equation_solver {
protected:
  func2var funct;
  double a, b;
  double cauchyProblem;
  uniform interval;
public:
  static double deviation( distribution const &valDistr1,
    distribution const &valDistr2 ) {

    return fabs(valDistr1[valDistr1.getNodeCount() - 1] - valDistr2[valDistr2.getNodeCount() - 1]);
  }

  virtual tabulated_function solve( double tollerance ) = 0;

  void setCauchyProblem( double cp ) { cauchyProblem = cp; }

  void setFragmentation( uint frag ) {
      interval.setBorders(frag + 1);
      interval.eval();
  }

  void setFunction( func2var f ) { funct = f; }

  void setBorders( double a, double b ) {
      this->a = a;
      this->b = b;
      interval.setBorders(a, b);
      interval.eval();
  }
};

class euler_cauchy_solver : public diff_equation_solver {
private: 
  uint minFrag;
  uint maxFrag;
public:
  euler_cauchy_solver( void ) {}

  static void eval( func2var funct, tabulated_function &s, double x1, double cp, double x2, uint localFrag ) {
    uniform argDistr(x1, x2);
    argDistr.setBorders(localFrag);
    argDistr.eval();
    double h = argDistr.getStep();
    s.clear();

    s << make_pair<>(x1, cp);

    for (uint i = 1; i < localFrag; i++) {
      auto y_pred = s[i - 1] + h * funct(argDistr[i - 1], s[i - 1]);
      s << make_pair<>(
               argDistr[i],
               s[i - 1] + (funct(argDistr[i - 1], s[i - 1]) + funct(argDistr[i], y_pred)) / 2 * h);
    }
  }

  tabulated_function solve( double tollerance ) override {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;

    solution << make_pair<>(a, cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + 1);

        for (; deviation(s2, s1) >= 1 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
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

class explicit_adams_solver : public diff_equation_solver {
public:
  explicit_adams_solver( void ) {}

  uint getMinFrag( void ) const { return minFrag; }
  uint getMaxFrag( void ) const { return maxFrag; }

  void eval( func2var funct, tabulated_function &s, double x1, double cp, double x2, uint localFrag ) {
    uniform argDistr(x1, x2);
    argDistr.setBorders(localFrag);
    argDistr.eval();
    double h = argDistr.getStep();
    s.clear();

    // Preset using RK2
    euler_cauchy_solver::eval(funct, s, x1, cp, x1 + h, 2);

    for (uint i = 2; i < localFrag; i++) {
      s << make_pair<>(
               argDistr[i],
               s[i - 1] + h / 2 * (3 * funct(argDistr[i - 1], s[i - 1]) - funct(argDistr[i - 2], s[i - 2])));
    }
  }

  tabulated_function solve( double tollerance ) override {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;

    solution << make_pair<>(a, cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + 1);

        for (; deviation(s2, s1) >= 3 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }

private:
  uint minFrag;
  uint maxFrag;
};

class implicit_adams_solver : public diff_equation_solver {
public:
  implicit_adams_solver( void ) {}

  uint getMinFrag( void ) const { return minFrag; }
  uint getMaxFrag( void ) const { return maxFrag; }

  void eval( func2var funct, tabulated_function &s, double x1, double cp, double x2, uint localFrag ) {
    uniform argDistr(x1, x2);
    argDistr.setBorders(localFrag);
    argDistr.eval();
    double h = argDistr.getStep();
    s.clear();

    // Preset using RK2
    euler_cauchy_solver::eval(funct, s, x1, cp, x1 + h, 2);

    for (uint i = 2; i <= localFrag; i++) {
      auto
        predictor = s[i - 1] + h / 2 * (3 * funct(argDistr[i - 1], s[i - 1]) - funct(argDistr[i - 2], s[i - 2])),
        corrector = s[i - 1] + h / 2 * (funct(argDistr[i], predictor) + funct(argDistr[i - 1], s[i - 1]));
      s << make_pair<>(argDistr[i], corrector);
    }
  }

  tabulated_function solve( double tollerance ) override {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;

    solution << make_pair<>(a, cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + 1);

        for (; deviation(s2, s1) >= 3 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }

private:
  uint minFrag;
  uint maxFrag;
};
