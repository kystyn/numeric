#include <iostream>
#include "controller.h"

using namespace std;

int main()
{
    controller c("integral.in");
    c << controller::nonDiffFunc << controller::diffFunc;
    c.run("integral.out");
    return 0;
}
