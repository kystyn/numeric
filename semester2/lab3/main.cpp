#include <iostream>
#include "controller.h"

using namespace std;

int main()
{
    rado_integral_n2 i;
    i.setFunction(controller::diffFunc);
    i.setBorders(0, 2 * 3.14159265358979323846);
    i.setTollerance(1e-8);
    std::cout << i() << std::endl;
    //controller c("integral.in");
    //c << controller::nonDiffFunc << controller::diffFunc;
    //c.run("integral.out");
    return 0;
}
