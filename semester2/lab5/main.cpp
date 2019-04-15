#include <iostream>
#include "controller.h"

using namespace std;

std::ostream & operator<<( std::ostream &os, tabulated_function const &tf ) {
  for (auto &v : tf.Y)
    os << v << ' ';
  os << endl;

  return os;
}

int main()
{
    euler_cauchy_solver s;
    auto tf = s.setFunction([]( double x, double y ) { return 2 * x; }).
      setBorders(0, 1).
      setCauchyProblem(0).
      solve(1e-3);
    return 0;
}
