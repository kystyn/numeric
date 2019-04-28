#include <iostream>
#include "controller.h"

using namespace std;

std::ostream & operator<<( std::ostream &os, tabulated_function const &tf ) {
  for (auto &c : tf.Coordinates)
    os << c.first << ' ' << c.second << endl;

  return os;
}

int main()
{
    controller c("de.in");
    c << []( double x, double y ) -> double { return /*exp(x) / x +*/ y; };
    c.run("de.out");

    return 0;
}
