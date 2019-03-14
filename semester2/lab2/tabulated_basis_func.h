#ifndef BASIS_FUNC_H
#define BASIS_FUNC_H

#include <vector>
#include <functional>

#include "def.h"

// ATTENTION: Weights array which is given in constructor should exist until destructor is called!
class tabulated_basis_func {
private:
    std::vector<double> Func;
    std::vector<double> const &Weights;

public:
    tabulated_basis_func( std::function<double(double, int)> &F, std::vector<double> const &Grid, uint BasisFuncNum,
                std::vector<double> const &Weights ) :
        Weights(Weights), Func(Grid.size()) {
        for (uint i = 0, size = Grid.size(); i < size; i++)
            Func[i] = F(Grid[i], BasisFuncNum);
    }

    tabulated_basis_func( std::vector<double> const &Func, std::vector<double> const &Weights ) :
        Func(Func), Weights(Weights) {}

    double operator*( tabulated_basis_func const &F ) {
        double res = 0;

        for (uint i = 0, dim = Weights.size(); i < dim; i++)
            res += Func[i] * Weights[i] * F[i];

        return res;
    }

    double operator[]( uint i ) const {
        return Func[i];
    }

    double & operator[]( uint i ) {
        return Func[i];
    }
};

#endif // BASIS_FUNC_H
