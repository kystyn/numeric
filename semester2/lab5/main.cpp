#include <iostream>
#include "controller.h"


using namespace std;

ostream & operator<<( ostream &os, vector<double> const &f ) {
  for (auto x : f)
    os << x << ' ';
  os << endl;

  return os;
}

std::ostream & operator<<( std::ostream &os, tabulated_function const &tf ) {
  for (auto &c : tf.Coordinates)
    os << c.first << ' ' << c.second << endl;

  return os;
}

vector<double> operator*( double h, vector<double> const &f ) {
  auto v = f;
  for (auto &x : v)
    x *= h;
  return v;
}

vector<double> operator*( vector<double> const &v1, vector<double> const &v2 ) {
  auto v = v1;
  for (uint i = 0; i < v2.size(); i++)
    v[i] *= v2[i];
  return v;
}

vector<double> operator-( vector<double> const &v1, vector<double> const &v2 ) {
  auto v = v1;

  for (uint i = 0; i < v2.size(); i++)
    v[i] -= v2[i];

  return v;
}

vector<double> operator+( vector<double> const &v1, vector<double> const &v2 ) {
  auto v = v1;

  for (uint i = 0; i < v2.size(); i++)
    v[i] += v2[i];
  return v;
}

vector<pair<double, vector<double>>> operator+( vector<pair<double, vector<double>>> const &v1,
                                                vector<pair<double, vector<double>>> const &v2 ) {
    auto v = v1;

    for (uint i = 0; i < v2.size(); i++)
      v[i].second = v[i].second + v2[i].second;
    return v;
}

vector<pair<double, vector<double>>> operator*( double h,
                                                vector<pair<double, vector<double>>> const &v ) {
    auto u = v;

    for (uint i = 0; i < v.size(); i++)
      u[i].second = h * v[i].second;
    return u;
}

// Chebyshev norm
double operator!( vector<double> const &v ) {
  double norm = 0;

  for (auto &x : v)
    if (fabs(x) > norm)
      norm = fabs(x);

  return norm;
}

int main()
{
    boundary_controller c("de.in");
    std::cout << std::setprecision(16) << kystyn::exp(1, 1e-14);

    c << array<func, 3>{[]( double x ) { return (2 * x + 2) / (2 * x * x + x);},
      []( double x ) { return -1.0 / (2 * x * x + x);},
      []( double x ) { return 1 / (x * (2 * x * x + x));}};
    c.run("de.out");

    /*
    reductor s;

    s.setBorders(uniform(0.2, 1, 5));
    s.setBoundaryProblem(1, 0, 5, 1, 0, 1);
    s.setFunction(
      []( double x ) { return (2 * x + 2) / (2 * x * x + x);},
      []( double x ) { return -1.0 / (2 * x * x + x);},
      []( double x ) { return 1 / (x * (2 * x * x + x));});
    auto sol = s.solve(1e-3);*/

    return 0;
}
