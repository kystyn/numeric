#pragma once

#include "matr.h"

namespace mth {
namespace lieqsys {
vec Relax( matr const &A, vec const &b, vec const &x0, double Omega, double Epsilon, int &Steps );
}
}
