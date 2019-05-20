#pragma once

#include <algorithm>

#include "def.h"
#include "distribution.h"
#include "interpolator.h"

class diff_equation_solver {
protected:
  uint minFrag;
  uint maxFrag;
  uint frag;
  uniform interval;
public:
  uint evalVolume;

  uint getMinFrag( void ) const { return minFrag; }
  uint getMaxFrag( void ) const { return maxFrag; }
  uint getFrag( void ) const { return frag; }

  static double deviation( tabulated_function const &valDistr1,
      tabulated_function const &valDistr2 ) {

      auto m = std::min(valDistr1[valDistr1.getNodeCount() - 1], valDistr2[valDistr2.getNodeCount() - 1]);

      //if (!m > 1e-14)
      //  return !((valDistr1[valDistr1.getNodeCount() - 1] - valDistr2[valDistr2.getNodeCount() - 1]) * m);
      return !(valDistr1[valDistr1.getNodeCount() - 1] - valDistr2[valDistr2.getNodeCount() - 1]);
  }

  virtual tabulated_function solve( double tollerance ) = 0;

  void setFragmentation( uint frag ) {
      interval.setBorders(frag + 1);
      interval.eval();
  }

  void setBorders( double a, double b ) {
      interval.setBorders(a, b);
      interval.eval();
  }

  void setBorders( uniform const &i ) {
    interval = i;
  }
};

class cauchy_problem_solver : public diff_equation_solver {
protected:
  n_func_n_plus_1_var funct;
  vector<double> cauchyProblem;

public:
  void setCauchyProblem( vector<double> cp ) { cauchyProblem = cp; }

  void setFunction( n_func_n_plus_1_var f ) { funct = f; }
};

class boundary_problem_solver : public diff_equation_solver {
protected:
  func p, q, f;
  double
    alpha0, alpha1, A,
    beta0, beta1, B;

public:
  void setBoundaryProblem( double a0, double a1, double A, double b0, double b1, double B ) { 
    alpha0 = a0;
    alpha1 = a1;
    this->A = A;

    beta0 = b0;
    beta1 = b1; 
    this->B = B;
  }

  void setFunction( func p1, func q1, func f1 ) { p = p1; q = q1; f = f1; }
};

class euler_cauchy_solver : public cauchy_problem_solver {
public:
  euler_cauchy_solver( void ) {}

  static void eval( n_func_n_plus_1_var funct, tabulated_function &s, double x1, vector<double> cp, double x2, uint localFrag ) {
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
               s[i - 1] + h / 2 * (funct(argDistr[i - 1], s[i - 1]) + funct(argDistr[i], y_pred)));
    }
  }

  tabulated_function solve( double tollerance ) override {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;
    frag = 0;
    evalVolume = 0;

    solution << make_pair<>(interval[0], cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + vector<double>{10 * tollerance});

        for (; deviation(s2, s1) >= 3 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
            evalVolume += 1 + (localFrag - 1) * 3;
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        frag += s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }
};

class explicit_adams_solver : public cauchy_problem_solver {
public:
  explicit_adams_solver( void ) {}

  static void eval( n_func_n_plus_1_var funct, tabulated_function &s, double x1, vector<double> cp, double x2, uint localFrag ) {
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
    frag = 0;

    solution << make_pair<>(interval[0], cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + vector<double>{10 * tollerance});

        for (; deviation(s2, s1) >= 3 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        frag += s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }
};

class implicit_adams_solver : public cauchy_problem_solver {
public:
  implicit_adams_solver( void ) {}

  void eval( n_func_n_plus_1_var funct, tabulated_function &s, double x1, vector<double> cp, double x2, uint localFrag ) {
    uniform argDistr(x1, x2);
    argDistr.setBorders(localFrag);
    argDistr.eval();
    double h = argDistr.getStep();
    s.clear();

    s << make_pair<>(x1, cp);

    for (uint i = 1; i < localFrag; i++) {
      tabulated_function RKsol;
      euler_cauchy_solver::eval(funct, RKsol, argDistr[i - 1], s[i - 1], argDistr[i], 3);
      auto yi = RKsol[2];
      s << make_pair<>(argDistr[i],
        s[i - 1] + h / 2 * (funct(argDistr[i], yi) + funct(argDistr[i - 1], s[i - 1])));

    }
  }

  tabulated_function solve( double tollerance ) override {
    tabulated_function solution;
    minFrag = 0xFFFFFFFF;
    maxFrag = 0;
    frag = 0;

    solution << make_pair<>(interval[0], cauchyProblem);

    for (uint i = 0; i < interval.getNodeCount() - 1; i++) {
        tabulated_function s1, s2;
        uint localFrag = 2;
        s2 << solution.get(solution.getNodeCount() - 1);
        auto p = s2.get(0);
        s1 << make_pair<>(p.first, p.second + vector<double>{10 * tollerance});

        for (; deviation(s2, s1) >= 3 * tollerance;) {
            s1 = s2;
            eval(funct, s2, interval[i], solution[i], interval[i + 1], localFrag = ((localFrag - 1) << 1) + 1);
        }

        if (s2.getNodeCount() < minFrag)
            minFrag = s2.getNodeCount();
        if (s2.getNodeCount() > maxFrag)
            maxFrag = s2.getNodeCount();

        frag += s2.getNodeCount();

        solution << s2.get(s2.getNodeCount() - 1);
    }
    return solution;
  }
};

class finite_difference_solver : public boundary_problem_solver {
public:
  static vector<double> tridiagSolve( vector<double> const &a, vector<double> const &b, vector<double> const &c, vector<double> const &d ) {
    vector<double> alpha(b.size()), beta(b.size()), x(b.size());
    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (uint i = 1; i < b.size(); i++) {
      alpha[i] = -c[i] / (a[i - 1] * alpha[i - 1] + b[i]);
      beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / (a[i - 1] * alpha[i - 1] + b[i]);
    }
    x[x.size() - 1] = beta[x.size() - 1];
    for (int i = x.size() - 2; i >= 0; i--)
      x[i] = alpha[i] * x[i + 1] + beta[i];

    return x;
  }

  static double tridiagSolveLast( vector<double> const &a, vector<double> const &b, vector<double> const &c, vector<double> const &d ) {
    double alpha, beta, x;
    alpha = -c[0] / b[0];
    beta = d[0] / b[0];
    for (uint i = 1; i < b.size(); i++) {
      auto oldalpha = alpha;
      alpha = -c[i] / (a[i - 1] * alpha + b[i]);
      beta = (d[i] - a[i - 1] * beta) / (a[i - 1] * oldalpha + b[i]);
    }
    x = beta;

    return x;
  }

  vector<double> eval( uniform segment ) {
    vector<double>
      a(segment.getNodeCount() - 1),
      b(segment.getNodeCount()),
      c(segment.getNodeCount(), 0),
      d(segment.getNodeCount());

    double
      h = segment.getStep();

    b[0] = h * alpha0 - alpha1;
    c[0] = alpha1;
    d[0] = A * h;

    a[a.size() - 1] = -beta1;
    b[b.size() - 1] = h * beta0 + beta1;
    d[d.size() - 1] = B * h;

    for (uint i = 1; i < b.size() - 1; i++) {
      a[i - 1] = 1 - h / 2 * p(segment[i]);
      b[i] = -(2 - h * h * q(segment[i]));
      c[i] = 1 + h / 2 * p(segment[i]);
      d[i] = h * h * f(segment[i]);
    }

    return tridiagSolve(a, b, c, d);
    //return tridiagSolveLast(a, b, c, d);
  }

 
  tabulated_function solve( double tollerance ) override {
      tabulated_function solution;
      frag = 0;
      evalVolume = 0;
    
      vector<double> s1, s2;
      uint localFrag = 1;
      uniform
        g1(interval.getA(), interval.getB(), (interval.getNodeCount() - 1) * localFrag + 1),
        g2(interval.getA(), interval.getB(), (interval.getNodeCount() - 1) * (localFrag <<= 1) + 1);
      s2 = eval(g1);

      bool isContinue = true;
      for (; isContinue;) {
          s1 = s2;
          s2 = eval(g2);
          evalVolume += g2.getNodeCount() * 4;
          isContinue = false;
          for (uint i = 0, j = 0; i < s1.size() && j < s2.size(); i += (s1.size() - 1) / (interval.getNodeCount() - 1), j += (s2.size() - 1) / (interval.getNodeCount() - 1)) {
              if (fabs(s1[i] - s2[j]) >= 3 * tollerance) {
                isContinue = true;
                break;
              }
          }
          localFrag <<= 1;
          g2.setBorders((interval.getNodeCount() - 1) * localFrag + 1);
          g2.eval();
      }

      frag = g2.getNodeCount();

      for (uint i = 0; i < s2.size(); i += (s2.size() - 1) / (interval.getNodeCount() - 1))
        solution << make_pair<>(interval[i / ((s2.size() - 1) / (interval.getNodeCount() - 1))], vector<double>{s2[i]});

      return solution;
}
};

class reductor : public boundary_problem_solver {
public:
  reductor( void ) {}

  tabulated_function solve( double tollerance ) override {
    euler_cauchy_solver ecs;
    ecs.setBorders(interval);
    ecs.setCauchyProblem({alpha1, -alpha0});
    ecs.setFunction([this]( double x, vector<double> const &y ) {
      std::vector<double> r(2);
      r[0] = y[1];
      r[1] = -p(x) * y[1] - q(x) * y[0];

      return r;
    });
    auto u = ecs.solve(tollerance);

    if (fabs(alpha0) >= tollerance)
        ecs.setCauchyProblem({A / alpha0, 0});
    else
        ecs.setCauchyProblem({0, A / alpha1});

    frag = ecs.getFrag();
    evalVolume = ecs.evalVolume;

    ecs.setFunction([this]( double x, vector<double> const &y ) {
      std::vector<double> r(2);
      r[0] = y[1];
      r[1] = f(x) - p(x) * y[1] - q(x) * y[0];

      return r;
    });
    auto v = ecs.solve(tollerance);

    frag = std::max(ecs.getFrag(), frag);
    evalVolume += ecs.evalVolume;
    minFrag = maxFrag = frag;

    double C = (B - beta0 * v[v.getNodeCount() - 1][0] - beta1 * v[v.getNodeCount() - 1][1]) /
            (beta0 * u[u.getNodeCount() - 1][0] + beta1 * u[u.getNodeCount() - 1][1]);

    auto solution = u * C + v;
    for (uint i = 0; i < solution.getNodeCount(); i++)
      solution[i].pop_back();


    return solution;
  }
};
