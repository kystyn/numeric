#include <iostream>
#include "def.h"
#include "controller.h"
#include "interpolator.h"

using namespace std;

int main()
{
    try {
        auto BasisFunc = [] ( double x, int j ) {
            //return cos(j * acos(0.6 + 0.5 * x));
            double r = 1;
            for (int i = 1; i <= j; i++)
                r *= x;
            return r;
        };

        controller c(BasisFunc, "lsm.in");

        auto normalDistribution = [] ( double mu, double sigma, double x ) {
            double e = (x - mu) / sigma;
            return exp(-e * e / 2) / sqrt(2 * pi * sigma);
        };

        c.
            addGrid('U', shared_ptr<uniform>(new uniform)).
            addGrid('R', shared_ptr<class random>(new class random)).
            addGrid('C', shared_ptr<chebyshev>(new chebyshev));


        c.addWeightDistribution('U', shared_ptr<constant>(new constant));

        auto normDistr = shared_ptr<normal>(new normal(normalDistribution));
        normDistr->makeBegin();
        c.addWeightDistribution('B', normDistr);

        normDistr = shared_ptr<normal>(new normal(normalDistribution));
        normDistr->makeMid();
        c.addWeightDistribution('M', normDistr);

        normDistr = shared_ptr<normal>(new normal(normalDistribution));
        normDistr->makeEnd();
        c.addWeightDistribution('E', normDistr);

        c << controller::nonDiffFunc << controller::diffFunc;


        c.run("lsm.out");
      } catch (const char *msg) {
        std::cerr << msg << std::endl;
      }

    return 0;
}
