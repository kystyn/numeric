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
    std::cout << std::setprecision(16) << kystyn::exp(1, 1e-14);

    c << []( double x, double y ) -> double { return kystyn::exp(x, 1e-14) / x + y; };
    c.run("de.out");

    return 0;
}
