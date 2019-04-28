#pragma once

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
struct data {
  double a, b;
  double cauchyProblem;
  uint fragmentation;
  double Tollerance;

  data( void ) : a(0), b(0), cauchyProblem(0),
      fragmentation(0), Tollerance(1e-16) {}

  data( double a, double b, double cauchyProblem, uint fragmentation,
        double Tollerance ) :
    a(a), b(b), cauchyProblem(cauchyProblem),
    fragmentation(fragmentation), Tollerance(Tollerance) {}
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
      uint frag;
      double tollerance;

      if (f >> a)
        f >> b >> cp >> frag >> tollerance;
      else
        break;
      loadedData.push_back(kyst::data(a, b, cp, frag, tollerance));
    }

    return *this;
  }

public:
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
      ecs.setBorders(d.a, d.b);
      ecs.setFunction(funct);
      ecs.setCauchyProblem(d.cauchyProblem);
      ecs.setFragmentation(d.fragmentation);
      auto tf = ecs.solve(d.Tollerance);
        fs << std::setprecision(16) << ecs.getMinFrag() << ' ' << ecs.getMaxFrag() << endl << tf << endl;
    }
  }
};
