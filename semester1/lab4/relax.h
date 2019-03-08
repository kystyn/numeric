#pragma once

#include "matr.h"

vec Relax( matr const &A, vec const &b, vec const &x0, double Omega, double Epsilon, int &Steps );
