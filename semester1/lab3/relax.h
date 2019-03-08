#pragma once

#include <vector>

std::vector<double> Relax( std::vector<std::vector<double>> const &A, std::vector<double> const &b,
                           std::vector<double> const &x0, double Omega, double Epsilon, int &Steps );
