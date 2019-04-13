#include <iostream>
#include "controller.h"

using namespace std;

std::ostream & operator<<( std::ostream &os, tabulated_function const &tf ) {
  for (auto &v : tf.Y)
    os << v << ' ';
  os << endl;

  return os;
}

int main()
{
    return 0;
}
