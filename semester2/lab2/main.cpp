#include <iostream>
#include "def.h"
#include "controller.h"

using namespace std;

int main()
{
    pow(1, 0);
    try {
        controller c("lsm.in", []( double x, int j ) { return pow(x, j); });


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
