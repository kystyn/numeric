#pragma once

#include <array>
#include <fstream>
#include <string>
#include <random>
#include <memory>
#include <vector>
#include <iomanip>

#include "distribution.h"
#include "diff_equations.h"

using namespace std;

namespace kyst {
struct boundary_data {
  uint methodN;
  double a, b;
  double alpha0, alpha1, A, beta0, beta1, B;
  uint fragmentation;
  double Tollerance;

  boundary_data( void ) : a(0), b(0), alpha0(0),
      fragmentation(0), Tollerance(1e-16) {}

  boundary_data( uint methodN, double a, double b, double alpha0, double alpha1, double A, double beta0, double beta1, double B, uint fragmentation,
        double Tollerance ) :
    methodN(methodN), a(a), b(b), alpha0(alpha0), alpha1(alpha1), beta0(beta0), beta1(beta1), A(A), B(B),
    fragmentation(fragmentation), Tollerance(Tollerance) {}
};

struct cauchy_data {
  uint methodN;
  double a, b;
  double cauchyProblem;
  uint fragmentation;
  double Tollerance;

  cauchy_data( void ) : a(0), b(0), cauchyProblem(0),
      fragmentation(0), Tollerance(1e-16) {}

  cauchy_data( uint methodN, double a, double b, double cauchyProblem, uint fragmentation,
        double Tollerance ) :
    methodN(methodN), a(a), b(b), cauchyProblem(cauchyProblem),
    fragmentation(fragmentation), Tollerance(Tollerance) {}
};
}

class cauchy_controller {
private:
  //euler_cauchy_solver ecs;
  std::shared_ptr<cauchy_problem_solver> cps;
  n_func_n_plus_1_var funct;
  vector<kyst::cauchy_data> loadedData;

  cauchy_controller & loadFromFile( const char *fileName ) {
    ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      uint methodN;
      double a, b;
      double cp;
      uint frag;
      double tollerance;

      if (f >> methodN)
        f >> a >> b >> cp >> frag >> tollerance;
      else
        break;
      loadedData.push_back(kyst::cauchy_data(methodN, a, b, cp, frag, tollerance));
    }

    return *this;
  }

public:
  cauchy_controller( const char *fileName = nullptr ) {
      if (fileName)
        loadFromFile(fileName);
  }

  cauchy_controller & operator<<( n_func_n_plus_1_var v ) {
    funct = v;
    return *this;
  }

  void run( const char *fileName ) {
    auto fs = ofstream(fileName);

    for (auto &d : loadedData) {
      switch (d.methodN) {
      case 0:
        cps = shared_ptr<euler_cauchy_solver>(new euler_cauchy_solver);
        break;
      case 1:
        cps = shared_ptr<explicit_adams_solver>(new explicit_adams_solver);
        break;
      case 2:
        cps = shared_ptr<implicit_adams_solver>(new implicit_adams_solver);
        break;
      }
      cps->setBorders(d.a, d.b);
      cps->setFunction(funct);
      cps->setCauchyProblem({d.cauchyProblem});
      cps->setFragmentation(d.fragmentation);
      auto tf = cps->solve(d.Tollerance);
        fs << std::setprecision(16) << cps->getFrag() << ' ' << cps->getMinFrag() << ' ' << cps->getMaxFrag() << endl << tf << endl;
    }
  }
};

class boundary_controller {
private:
  std::shared_ptr<boundary_problem_solver> bps;
  func p, q, f;
  vector<kyst::boundary_data> loadedData;

  boundary_controller & loadFromFile( const char *fileName ) {
    ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      uint methodN;
      double a, b;
      double alpha0, alpha1, beta0, beta1, A, B;
      uint frag;
      double tollerance;

      if (f >> methodN)
        f >> a >> b >> alpha0 >> alpha1 >> A >> beta0 >> beta1 >> B >> frag >> tollerance;
      else
        break;
      loadedData.push_back(kyst::boundary_data(methodN, a, b, alpha0, alpha1, A, beta0, beta1, B, frag, tollerance));
    }

    return *this;
  }

public:
  boundary_controller( const char *fileName = nullptr ) {
      if (fileName)
        loadFromFile(fileName);
  }

  boundary_controller & operator<<( array<func, 3> f ) {
    p = f[0];
    q = f[1];
    this->f = f[2];
    return *this;
  }

  void run( const char *fileName ) {
    auto fs = ofstream(fileName);

    for (auto &d : loadedData) {
      switch (d.methodN) {
      case 0:
        bps = shared_ptr<finite_difference_solver>(new finite_difference_solver);
        break;
      case 1:

        break;
      }
      bps->setBorders(d.a, d.b);
      bps->setFunction(p, q, f);
      bps->setBoundaryProblem(d.alpha0, d.alpha1, d.A, d.beta0, d.beta1, d.B);
      bps->setFragmentation(d.fragmentation);
      auto tf = bps->solve(d.Tollerance);
        fs << std::setprecision(16) << bps->getFrag() << ' ' << bps->getMinFrag() << ' ' << bps->getMaxFrag() << endl << tf << endl;
    }
  }
};
