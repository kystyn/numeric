#include <iostream>
#include "controller.h"

using namespace std;

int main()
{
    //trapezium_integral i;
    //i.setFunction([]( double x ) { return x; });
    //i.setBorders(0, 2);
    //i.setTollerance(1e-8);
    //std::cout << i() << std::endl;
    controller c("integral.in");
    c << controller::nonDiffFunc << controller::diffFunc;
    c.run("integral.out");
    return 0;
}
