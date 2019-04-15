#include <iostream>
#include "controller.h"

using namespace std;

std::ostream & operator<<( std::ostream &os, tabulated_function const &tf ) {
  for (uint i = 0; i < tf.Y.size(); i++)
    os << (*tf.X)[i] << ' ' << tf.Y[i] << endl;
  os << endl;

  return os;
}

int main()
{
   /*
  euler_cauchy_solver s;
    auto tf = s.setFunction([]( double x, double y ) { return exp(x) * (log(x) + 1); }).
      setBorders(1, 3).
      setCauchyProblem(exp(1)).
      solve(1e-4);
    */
  
    controller c("de.in");
    c << []( double x, double y ) -> double { return exp(x) / x + y; };
    c.run("de.out");
    
    return 0;
}
