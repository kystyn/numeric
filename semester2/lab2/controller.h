#pragma once

#include <map>
#include <fstream>
#include <string>
#include <random>
#include "distribution.h"
#include "interpolator.h"

struct data {
  double a, b;
  int nodeCount;
  char gridT;
  int weightsT;
  int weightsCount;
  int funcT;

  data( void ) : a(0), b(0), nodeCount(0),
      gridT('u'), weightsT(0), weightsCount(0), funcT(0){}

  data( double a, double b, int node, char gr, int weightsT, int weightsCnt, int f ) :
    a(a), b(b), nodeCount(node),
    gridT(gr), weightsT(weightsT),
    weightsCount(weightsCnt), funcT(f) {}
};

class controller {
private:
  least_square_interpolator i;

  std::map<char, distribution *> availableGridX;
  std::vector<func> availableFunctions;
  std::vector<data> loadedData;

  controller & loadFromFile( const char *fileName ) {
    std::ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      double a, b;
      char gridT;
      int weightsT, funcT, weightsCount;
      int nc;
      if (f >> a)
        f >> b >> nc >> gridT >> weightsT >> weightsCount >> funcT;
      else
        break;
      loadedData.push_back(data(a, b, nc, gridT, weightsT, weightsCount, funcT));
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
    return fabs(a) < mth::Tollerance ? 1 : -1;
  }

  static double derivativeOfNonDiffFunc( double x ) {
    return (double)pow_1(x < 3) * -sin(2 * x) / (1 + cos(x) * cos(x));
    //return 8 * pow(x, 7);
  }

  controller( least_square_interpolator::basis_func &F ) {
      i.setBasisFunc(F);
  }

  distribution *& operator[]( char p ) {
    return availableGridX[p];
  }

  controller & operator<<( func v ) {
    availableFunctions.push_back(v);
    return *this;
  }

  controller( const char *fileName, least_square_interpolator::basis_func F ) {
    loadFromFile(fileName);
    i.setBasisFunc(F);
  }

  void run( const char *fileName ) {
    auto fs = std::ofstream(fileName);

    for (auto &d : loadedData) {
      auto f = availableFunctions[d.funcT];
      i.setFunc(f);

      std::vector<double> Weights(d.nodeCount);

      switch (d.weightsT) {
      case 0:
          Weights = std::vector<double>(d.nodeCount, 1);
          break;
      case 1: {
          std::default_random_engine RandomEngine;
          std::normal_distribution<double> distr(1, 3);
          for (auto &w : Weights)
              w = distr(RandomEngine);
          }
          break;
      case 2: {
          std::default_random_engine RandomEngine;
          std::exponential_distribution<double> distr(2);
          for (auto &w : Weights)
              w = 1 + distr(RandomEngine);
          }
          break;
      }
      i.setWeights(Weights).setGrid(availableGridX[d.gridT]->setBorders(d.a, d.b, d.nodeCount).eval()).genPolinom();
      fs << d.gridT << ' ' << d.funcT << ' ' << d.weightsT << '\n';
      i.output(fs);
    }
  }
};
