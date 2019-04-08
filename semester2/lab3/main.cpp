#include <iostream>
#include "controller.h"

using namespace std;

int main()
{
    rado_integral_n2 i;
    i.setFunction(controller::diffFunc);
    i.setGrid(uniform(0, 1.5));
    i.setTollerance(1e-6);
    std::cout << i() << std::endl;
    //controller c("integral.in");
    //c << controller::nonDiffFunc << controller::diffFunc;
    //c.run("integral.out");
    return 0;
}
