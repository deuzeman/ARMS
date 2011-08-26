#include <iostream>
#include <Data.h>
#include <Minim.h>

int main(int argc, char **argv)
{
  Data data("test.dat");
  Minim minim(data, 20, 1, 12, 234231);
  Eigen::ArrayXd start(4);
  start << 0.05, 0.0, 0.0, 0.02;
  Eigen::ArrayXXd bounds(2,4);
  bounds << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0; 
  minim.powell(start, bounds, 50, 5, 1e-8);
  return 0;
}
