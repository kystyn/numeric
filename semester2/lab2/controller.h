#pragma once

#include <map>
#include <fstream>
#include <string>
#include <random>
#include <memory>
#include "distribution.h"
#include "interpolator.h"

struct data {
  double a, b;
  uint nodeCount;
  char gridT;
  char weightsT;
  uint funcT;
  uint polynomDegree;

  data( void ) : a(0), b(0), nodeCount(0),
      gridT('u'), weightsT(0), funcT(0), polynomDegree(0) {}

  data( double a, double b, uint node, char grT, char weightsT, uint funcT, uint polynomDegree ) :
    a(a), b(b), nodeCount(node),
    gridT(grT), weightsT(weightsT),
    funcT(funcT), polynomDegree(polynomDegree) {}
};

class controller {
private:
  least_square_interpolator i;

  map<char, shared_ptr<distribution>> availableGridX;
  map<char, shared_ptr<distribution>> availableWeights;
  vector<func> availableFunctions;
  vector<data> loadedData;

  controller & loadFromFile( const char *fileName ) {
    ifstream f(fileName);
  
    for (int i = 0; ; i++) {
      double a, b;
      char gridT, weightsT;
      uint funcT;
      uint nodeCount;
      uint polynomDegree;

      if (f >> a)
        f >> b >> nodeCount >> gridT >> weightsT >> funcT >> polynomDegree;
      else
        break;
      loadedData.push_back(data(a, b, nodeCount, gridT, weightsT, funcT, polynomDegree));
    }

    return *this;
  }

public:

  static double diffFunc( double x ) {
    return log(1 + cos(x) * cos(x));
    //return pow(x, 8) + 1;
  }

  static double nonDiffFunc( double x ) {
    return fabs(log(1 + cos(x) * cos(x)) - 0.4);
    //return pow(x, 8) + 1;
  }

  static double pow_1( double a ) {
    return fabs(a) < mth::Tollerance ? 1 : -1;
  }

  controller( least_square_interpolator::basis_func F, const char *fileName = nullptr ) {
      i.setBasisFunc(F);

      if (fileName)
        loadFromFile(fileName);
  }

  controller & addGrid( char p, shared_ptr<distribution> const distr ) {
      availableGridX[p] = distr;
      return *this;
  }

  controller & addWeightDistribution( char p, shared_ptr<distribution> const distr ) {
      availableWeights[p] = distr;
      return *this;
  }


  controller & operator<<( func v ) {
    availableFunctions.push_back(v);
    return *this;
  }

  void run( const char *fileName ) {
    auto fs = ofstream(fileName);

    for (auto &d : loadedData) {
        auto f = availableFunctions[d.funcT];
        i.
            setPolyDimension(d.polynomDegree).
            setFunc(f).
            setGrid(availableGridX[d.gridT]->setBorders(d.a, d.b, d.nodeCount).eval()).
            setWeights(availableWeights[d.weightsT]->setBorders(d.a, d.b, d.nodeCount).eval()).
            genPolinom();

      i.output(fs);
    }
  }
};
