#ifndef DEF_H
#define DEF_H

#include <functional>

using namespace std;

using func = function<double(double)>;
using func2var = function<double(double, double)>;

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

#endif // DEF_H
