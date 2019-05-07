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

  virtual distribution & setBorders( double newA, double newB ) {
    a = newA;
    b = newB;

    return *this;
  }

  virtual distribution & setBorders( size_t newNodeCount ) {
      NodeCount = newNodeCount;

      return *this;
  }

  virtual double operator[]( uint i ) const { return i; }

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

  double operator[]( int i ) const {
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

class grid : public distribution {
public:
  grid( std::vector<double> const &v = {} ) : Values(v) {
    NodeCount = v.size();
  }

  grid & operator<<( double val ) {
    Values.push_back(val);
    a = Values[0];
    b = val;
    NodeCount = Values.size();
    return *this;
  }

  double operator[]( uint i ) const {
    return Values[i];
  }
private:
  vector<double> Values;
};

class value_distribution : public distribution {
private:
  func F;
  shared_ptr<distribution> DistrX;
public:
  value_distribution( void ) : F(nullptr), DistrX(nullptr) {}

  value_distribution( shared_ptr<distribution> DistrX, func F ) :
    distribution(DistrX->getA(), DistrX->getB(), DistrX->getNodeCount()), F(F), DistrX(DistrX) { /*f (F != nullptr) eval();*/ }

  value_distribution & setGridX( shared_ptr<distribution> NewDistrX ) {
      DistrX = NewDistrX;
      return *this;
  }

  shared_ptr<distribution> X() { return DistrX; }

  value_distribution & setFunc( func Fun ) {
      F = Fun;
      return *this;
  }

  distribution & eval( void ) {
      DistrX->eval();
      return *this;
  }

  double operator[]( uint i ) const {
      if (F == nullptr)
          throw "Set function!";
      return F((*DistrX)[i]);
  }

  double operator()( double X ) const {
      return F(X);
  }

  ~value_distribution( void ) {}
};

class tabulated_function {
private:
    vector<pair<double, vector<double>>> Coordinates;
    uint NodeCount;
public:
    tabulated_function( void ) {}

    tabulated_function( vector<pair<double, vector<double>>> Coordinates ) : Coordinates(Coordinates) {}

    uint getNodeCount( void ) const { return NodeCount; }

    tabulated_function & operator<<( pair<double, vector<double>> p ) {
      Coordinates.push_back(p);
      NodeCount = Coordinates.size();
      return *this;
    }

    vector<double> operator[]( uint idx ) const {
      return Coordinates.at(idx).second;
    }

    pair<double, vector<double>> get( uint idx ) const {
        return Coordinates.at(idx);
    }

    tabulated_function & clear( void ) {
      Coordinates.clear();
      NodeCount = 0;
      return *this;
    }

    tabulated_function & pop( void ) {
      Coordinates.pop_back();
      NodeCount--;
      return *this;
    }

    tabulated_function & reverse( void ) {
      std::reverse(Coordinates.begin(), Coordinates.end());
      return *this;
    }

    tabulated_function operator+( tabulated_function const &tf ) const {
        return tabulated_function(Coordinates + tf.Coordinates);
    }

    tabulated_function operator*( double h ) const {
        return tabulated_function(h * Coordinates);
    }

    friend std::ostream & operator<<( std::ostream &os, tabulated_function const &tf );
};
