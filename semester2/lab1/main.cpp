#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include "interpolator.h"
#include "distribution.h"
#include "controller.h"


int main( void ) {

  try {
    //const double pi = acos(-1);
    controller c("hermit.in");
    std::vector<double (*)(double)> v;

    c['U'] = new uniform();
    c['R'] = new class random();
    c['C'] = new chebyshev();

    c << controller::nonDiffFunc << controller::diffFunc;

    v.push_back(controller::derivativeOfNonDiffFunc);
    c << v;

    v[0] = controller::derivativeOfDiffFunc;
    c << v;

    c.run("hermit.out");
//    interpolator i;
//    i.setFunc(sin, {cos});
//    i.setGrid(distribution(std::vector<double>{0.0, pi / 6, pi / 2})).buildDividedDifferenceTable();
  } catch (const char *msg) {
    std::cerr << msg << std::endl;
  }

  return 0;
}
