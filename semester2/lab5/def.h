#ifndef DEF_H
#define DEF_H

#include <functional>
#include <vector>


class tabulated_function;

using namespace std;

using func = function<double(double)>;
using n_func_n_plus_1_var = function<vector<double>(double, vector<double> const &)>;

using uint = unsigned int;

static const double pi = 3.14159265358979323;

namespace kystyn {
inline double exp( double x, double tollerance ) {
  long double cur = 1;
  long double res = cur;
  uint n = 1;

  while (fabs(cur) > tollerance) {
    cur *= x / n;
    n++;
    res += cur;
  }
  return res;
}

inline double log( double x, double tollerance ) {
  x = x - 1;
  long double cur = x;
  long double res = cur;
  uint n = 1;

  while (fabs(cur) > tollerance) {
    cur *= -x / n;
    n++;
    res += cur;
  }
  return res;
}
}

ostream & operator<<( ostream &os, vector<double> const &f );
std::ostream & operator<<( std::ostream &os, tabulated_function const &tf );
vector<double> operator*( double h, vector<double> const &f );
vector<double> operator*( vector<double> const &v1, vector<double> const &v2 );
vector<double> operator-( vector<double> const &v1, vector<double> const &v2 );
vector<double> operator+( vector<double> const &v1, vector<double> const &v2 );
vector<pair<double, vector<double>>> operator+( vector<pair<double, vector<double>>> const &v1,
                                                vector<pair<double, vector<double>>> const &v2 );
vector<pair<double, vector<double>>> operator*( double h,
                                                vector<pair<double, vector<double>>> const &v );

// Chebyshev norm
double operator!( vector<double> const &v );

#endif // DEF_H
