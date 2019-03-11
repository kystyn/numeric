#pragma once

#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include <cmath>

#include "def.h"

class interpolator;

class distribution {
public:
  double a, b;
  size_t NodeCount;
  bool NeedEvalGrid;
  std::vector<double> Grid;

  distribution( void ) {}

  distribution( std::vector<double> const &gr ) : a(gr.front()), b(gr.back()), NodeCount(gr.size()), Grid(gr) {}

  distribution( double a, double b, size_t NodeCount ) : a(a), b(b), NodeCount(NodeCount), NeedEvalGrid(true) { Grid.resize(NodeCount); }
  std::vector<double> const & operator()( void ) { return Grid; }

  virtual distribution & eval( void ) { return *this; }

  distribution & setBorders( double newA, double newB, size_t newNodeCount ) {
    a = newA;
    b = newB;
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
  uniform( double a, double b, size_t NodeCount ) : distribution(a, b, NodeCount) {}

  distribution & eval( void ) {
    double h = (b - a) / (double)(NodeCount - 1);
    int step = 0;
    for (std::vector<double>::iterator it = Grid.begin(); it < Grid.end(); ++it)
      *it = a + step++ * h;

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
      pi = acos(-1);
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
    distribution(DistrX.a, DistrX.b, DistrX.NodeCount), F(F), DistrX(DistrX) {}

  distribution & eval( void ) {
    std::vector<double>::const_iterator itX;
    std::vector<double>::iterator itY;
    for (itX = DistrX().begin(),
         itY = Grid.begin(); itX < DistrX().end(); ++itX, ++itY)
      *itY = F(*itX);

    return *this;
  }

  ~value_distribution( void ) {}
};
