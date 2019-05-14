#include <iostream>
#include "controller.h"

using namespace std;

int main()
{
    //rado_integral_n2 i;
    //i.setFunction([]( double x ) { return x * x * x * x; });
    //i.setBorders(0, 2);
    //i.setTollerance(1e-8);
    //std::cout << i() << std::endl;
    controller c("integral.in");
    c << controller::nonDiffFunc << controller::diffFunc;
    c.run("integral.out");
    return 0;
}
