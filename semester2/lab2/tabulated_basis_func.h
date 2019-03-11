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
    tabulated_basis_func( std::function<double(double, int)> &F, double ValX,
                std::vector<double> const &Weights ) :
        Weights(Weights), Func(Weights.size()) {
        for (uint i = 0, size = Weights.size(); i < size; i++)
            Func[i] = F(ValX, i);
    }

    tabulated_basis_func( std::vector<double> const &F, std::vector<double> const &W ) : Func(F), Weights(W) {
        if (F.size() != W.size())
            throw "Bad basis function!";
    }

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
