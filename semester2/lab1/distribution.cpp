#include "distribution.h"
#include "interpolator.h"

/*bad::bad( double a, double b, int nodeCount, interpolator &in ) : distribution(a, b, nodeCount), interpol(in) {}

distribution & bad::eval( void ) {
  const double precision = 1e2;
  distribution uniformGrid = uniform(a, b, (b - a) * precision);

  interpol.setGrid(uniformGrid).buildDividedDifferenceTable();

  std::vector<std::pair<double, double>> v(nodeCount);

  int i = 0;
  for (std::vector<std::pair<double, double>>::iterator it = v.begin(); it < v.end(); ++it, i++) {
    it->first = uniformGrid[i];
    it->second = fabs(interpol.HermitPolinom(it->first) - interpol.getFunc(it->first));
  }

  std::sort(v.begin(), v.end(), &compare);

  grid.resize(nodeCount);

  for (int i = 0; i < nodeCount; i++)
    grid[i] = v[i].first;

  std::sort(grid.begin(), grid.end());

  needEvalGrid = false;

  return *this;
}*/