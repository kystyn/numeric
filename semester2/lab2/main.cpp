#include <iostream>
#include "def.h"
#include "controller.h"
#include "interpolator.h"

using namespace std;

int main()
{
    pow(1, 0);
    try {
        constexpr int dim = 3;
        constexpr int nodes = 11;
        auto Func = [] ( double x ) {
            return sin(x * x * x);
        };
        auto BasisFunc = [] ( double x, int j ) {
            //return cos(j * acos(0.5 * (1 + x)));
            return pow(x, j);
        };
        least_square_interpolator i(Func, BasisFunc, dim, uniform(1, 2, nodes));
        i.setWeights(std::vector<double>(nodes, 1)).genPolinom();

        std::cout << "Deviation = " << !i << std::endl;

        controller c([]( double x, int j ) { return pow(x, j); }, dim, "lsm.in");


        c['U'] = new uniform();
        c['R'] = new class random();
        c['C'] = new chebyshev();

        c << controller::nonDiffFunc << controller::diffFunc;


        c.run("lsm.out");
      } catch (const char *msg) {
        std::cerr << msg << std::endl;
      }

    return 0;
}
