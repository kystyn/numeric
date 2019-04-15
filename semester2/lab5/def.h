#ifndef DEF_H
#define DEF_H

#include <functional>

using namespace std;

using func = function<double(double)>;
using func2var = function<double(double, double)>;

using uint = unsigned int;

static const double pi = 3.14159265358979323;

#endif // DEF_H
