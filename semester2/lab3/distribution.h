#pragma once

#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include <memory>
#include <cmath>

#include "def.h"

class distribution {
protected:
  double a, b;
  size_t NodeCount;

public:
  distribution( void ) : a(0), b(0), NodeCount(1) {}

  distribution( double a, double b, size_t NodeCount ) : a(a), b(b), NodeCount(NodeCount) {}

  virtual distribution & eval( void ) { return *this; }

  double getA( void ) const { return a; }
  double getB( void ) const { return b; }
  size_t getNodeCount( void ) const { return NodeCount; }

  distribution & setBorders( double newA, double newB ) {
    a = newA;
    b = newB;

    return *this;
  }

  distribution & setBorders( size_t newNodeCount ) {
      NodeCount = newNodeCount;

      return *this;
  }

  virtual double operator[]( uint i ) const { return 0; }


  virtual ~distribution() {}
};

class uniform : public distribution {
private:
    double Step;
public:
  uniform( void ) : Step(0) {}
  uniform( double a, double b, size_t NodeCount = 1 ) : distribution(a, b, NodeCount) {
      eval();
  }

  double operator[]( uint i ) const {
      return a + i * Step;
  }

  double getStep( void ) const {
      return Step;
  }

  distribution & eval( void ) {
    Step = (b - a) / (NodeCount - 1);
    return *this;
  }
};

class constant : public distribution {
public:
  constant( size_t NodeCount = 1 ) : distribution(1, 1, NodeCount) { eval(); }

  distribution & eval( void ) {
    return *this;
  }

  double operator[]( uint i ) const {
      return 1;
  }
};

class random : public distribution {
private:
    std::vector<double> Grid;
    std::default_random_engine RandomEngine;
    std::uniform_real_distribution<double> UniformDistr;
public:
  random( void ) {}
  random( double a, double b, size_t NodeCount ) : distribution(a, b, NodeCount), UniformDistr(0.5, 0.17) {}

  distribution & eval( void ) {
    double h = (b - a) / (double)(NodeCount - 1);
    int step = 0;
    for (std::vector<double>::iterator it = Grid.begin(); it < Grid.end(); ++it)
      *it = a + step++ * h;
    for (std::vector<double>::iterator it = Grid.begin(); it < Grid.end(); ++it)
      *it += h * UniformDistr(RandomEngine);

    return *this;
  }

  double operator[]( uint i ) const {
      return Grid[i];
  }
};

class chebyshev : public distribution {
public:
  chebyshev( void ) {}
  chebyshev( double a, double b, size_t NodeCount ) : distribution(a, b, NodeCount) {}

  distribution & eval( void ) {
    return *this;
  }

  double operator[]( uint i ) const {
      double
        sh = (a + b) / 2.0,
        mul = (b - a) / 2.0;

      return sh + mul * cos(pi * (2 * i + 1) / (2 * (NodeCount + 1)));
  }

  ~chebyshev( void ) {}
};

class value_distribution : public distribution {
private:
  func F;
  distribution DistrX;
public:
  value_distribution( void ) : F(nullptr) {}

  value_distribution( distribution const &DistrX, func F ) :
    distribution(DistrX.getA(), DistrX.getB(), DistrX.getNodeCount()), F(F), DistrX(DistrX) { /*f (F != nullptr) eval();*/ }

  value_distribution & setGridX( distribution const &NewDistrX ) {
      DistrX = NewDistrX;
      return *this;
  }

  distribution & setGridX( void ) {
      return DistrX;
  }

  value_distribution & setFunc( func Fun ) {
      F = Fun;
      return *this;
  }

  distribution & eval( void ) {
      return *this;
  }

  double operator[]( uint i ) const {
      if (F == nullptr)
          throw "Set function!";
      return F(DistrX[i]);
  }

  double operator()( double X ) const {
      return F(X);
  }

  ~value_distribution( void ) {}
};


