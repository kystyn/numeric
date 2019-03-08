#pragma once

#include <vector>
#include <ctime>
#include <algorithm>


class interpolator;

class distribution {
public:
  double a, b;
  int nodeCount;
  int needEvalGrid;
  std::vector<double> grid;

  distribution( void ) {}
  distribution( std::vector<double> const &gr ) : grid(gr), a(gr.front()), b(gr.back()), nodeCount(gr.size()) {}
  distribution( double a, double b, int nodeCount ) : a(a), b(b), nodeCount(nodeCount), needEvalGrid(true) { grid.resize(nodeCount); }
  std::vector<double> const & operator()( void ) { return grid; }

  virtual distribution & eval( void ) { return *this; }

  distribution & setBorders( double newA, double newB, int newNodeCount ) {
    a = newA;
    b = newB;
    nodeCount = newNodeCount;
    grid.resize(nodeCount);

    return *this;
  }

  double operator[]( int i ) const {
    return grid[i];
  }
};

class uniform : public distribution {
public:
  uniform( void ) {}
  uniform( double a, double b, int nodeCount ) : distribution(a, b, nodeCount) {}

  distribution & eval( void ) {
    double h = (b - a) / (double)(nodeCount - 1);
    int step = 0;
    for (std::vector<double>::iterator it = grid.begin(); it < grid.end(); ++it)
      *it = a + step++ * h;

    needEvalGrid = false;

    return *this;
  }
};

class random : public distribution {
public:
  random( void ) {}
  random( double a, double b, int nodeCount ) : distribution(a, b, nodeCount) {}

  distribution & eval( void ) {
    srand(time(NULL));
    double h = (b - a) / (double)(nodeCount - 1);
    int step = 0;
    for (std::vector<double>::iterator it = grid.begin(); it < grid.end(); ++it)
      *it = a + step++ * h;
    for (std::vector<double>::iterator it = grid.begin(); it < grid.end(); ++it)
      *it += h / (2.0 * (1 + rand() % 4));
    needEvalGrid = false;

    return *this;
  }
};

class chebyshev : public distribution {
public:
  chebyshev( void ) {}
  chebyshev( double a, double b, int nodeCount ) : distribution(a, b, nodeCount) {}

  distribution & eval( void ) {
    srand(time(NULL));
    double
      pi = acos(-1);
    double
      sh = (a + b) / 2.0,
      mul = (b - a) / 2.0;
    double step = 0;
    for (std::vector<double>::iterator it = grid.begin(); it < grid.end(); ++it)
      *it = sh + mul * cos(pi * (2 * step++ + 1) / (2 * (nodeCount + 1)));
    needEvalGrid = false;

    return *this;
  }
};
/*
class bad : public distribution {
private:
  interpolator &interpol;

  static bool compare( std::pair<double, double> const &val1,
    std::pair<double, double> const &val2 ) {
    return val1.first < val2.first;
  }
public:
  bad( void ) {}
  bad( double a, double b, int nodeCount, interpolator &in );
  distribution & eval( void );
};*/

class value_distribution : public distribution {
private:
  double (*func)( double );
  distribution distrX;
public:
  value_distribution( void ) : func(NULL) {} 

  value_distribution( distribution const &distrX, double (*func)( double ) ) : 
    distribution(distrX.a, distrX.b, distrX.nodeCount), func(func), distrX(distrX) {}

  distribution & eval( void ) {
    std::vector<double>::const_iterator itX;
    std::vector<double>::iterator itY;
    for (itX = distrX().begin(),
         itY = grid.begin(); itX < distrX().end(); ++itX, ++itY)
      *itY = func(*itX);

    return *this;
  }
};
