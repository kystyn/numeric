#pragma once

#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include <memory>
#include <cmath>

#include "def.h"

class distribution {
public:
  double a, b;
  size_t NodeCount;
  bool NeedEvalGrid;
  std::vector<double> Grid;

  distribution( void ) {}

  distribution( std::vector<double> const &gr ) : a(gr.front()), b(gr.back()), NodeCount(gr.size()), Grid(gr) {}

  distribution( double a, double b, size_t NodeCount ) : a(a), b(b), NodeCount(NodeCount), NeedEvalGrid(true) {
      Grid.resize(NodeCount);
  }
  std::vector<double> const & operator()( void ) { return Grid; }

  virtual distribution & eval( void ) { return *this; }

  distribution & setBorders( double newA, double newB, size_t newNodeCount = 1 ) {
    a = newA;
    b = newB;
    NodeCount = newNodeCount;
    Grid.resize(NodeCount);

    NeedEvalGrid = true;

    return *this;
  }

  distribution & setBorders( size_t newNodeCount ) {
      NodeCount = newNodeCount;
      Grid.resize(NodeCount);

       return *this;
  }

  double operator[]( uint i ) const {
    return Grid[i];
  }

  operator std::vector<double> const &( void ) const {
      return Grid;
  }

  virtual ~distribution() {}
};

class uniform : public distribution {
public:
  uniform( void ) {}
  uniform( double a, double b, size_t NodeCount = 1 ) : distribution(a, b, NodeCount) { eval(); }

  distribution & eval( void ) {
    double h = (b - a) / (double)(NodeCount - 1);
    int step = 0;
    for (std::vector<double>::iterator it = Grid.begin(); it < Grid.end(); ++it)
      *it = a + step++ * h;

    NeedEvalGrid = false;

    return *this;
  }
};

class constant : public distribution {
public:
  constant( size_t NodeCount = 1 ) : distribution(1, 1, NodeCount) { eval(); }

  distribution & eval( void ) {
    for (auto &g : Grid)
      g = 1;

    NeedEvalGrid = false;

    return *this;
  }
};

class random : public distribution {
private:
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
    NeedEvalGrid = false;

    return *this;
  }
};

class chebyshev : public distribution {
public:
  chebyshev( void ) {}
  chebyshev( double a, double b, size_t NodeCount ) : distribution(a, b, NodeCount) {}

  distribution & eval( void ) {
    double
      sh = (a + b) / 2.0,
      mul = (b - a) / 2.0;
    double step = 0;
    for (std::vector<double>::iterator it = Grid.begin(); it < Grid.end(); ++it)
      *it = sh + mul * cos(pi * (2 * step++ + 1) / (2 * (NodeCount + 1)));
    NeedEvalGrid = false;

    return *this;
  }

  ~chebyshev( void ) {}
};

class value_distribution : public distribution {
private:
  func F;
  distribution DistrX;
public:
  value_distribution( void ) {}

  value_distribution( distribution const &DistrX, func F ) :
    distribution(DistrX.a, DistrX.b, DistrX.NodeCount), F(F), DistrX(DistrX) { eval(); }

  distribution & eval( void ) {

    if (F == nullptr)
        throw "Function should be initialized first";

    std::vector<double>::const_iterator itX;
    std::vector<double>::iterator itY;
    for (itX = DistrX().begin(),
         itY = Grid.begin(); itX < DistrX().end(); ++itX, ++itY)
      *itY = F(*itX);

    return *this;
  }

  ~value_distribution( void ) {}
};

class normal : public distribution {
private:
    function<double(double, double, double)> Function;
    double *mu;
    double sigma;
    double mid;
public:
  normal( void ) : mu(&mid), sigma(1) {}
  normal( function<double(double, double, double)> F ) : Function(F), mu(&a), sigma(1) {}
  normal( double a, double b, size_t NodeCount ) : distribution(a, b, NodeCount), mu(&this->a), sigma(1) { eval(); }

  normal & makeBegin( void ) {
      mu = &a;
      sigma = 1;
      return *this;
  }

  normal & makeMid( void ) {
      mu = &mid;
      mid = (a + b) / 2.0;
      sigma = 1;
      return *this;
  }

  normal & makeEnd( void ) {
      mu = &b;
      sigma = 1;
      return *this;
  }

  distribution & eval( void ) {

    double h = (b - a) / (double)(NodeCount - 1);
    int i = 1;

    mid = (a + b) / 2.0;

    for (auto &g : Grid)
      g = Function(*mu, sigma, h * i++);

    NeedEvalGrid = false;

    return *this;
  }
};

