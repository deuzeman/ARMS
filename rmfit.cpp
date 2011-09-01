#include <iostream>
#include <Data.h>
#include <Minim.h>

int main(int argc, char **argv)
{
  Data data("test.dat");
  Minim minim(data, 20, 1, data.minEv(), data.maxEv(), 234231);
  Eigen::ArrayXd start(2);
  start << 0.07, 0.05;
  Eigen::ArrayXXd bounds(2,2);
  bounds << 0.0, 0.0, 0.5, 0.15; 
  minim.powell(start, bounds, 5000, 50, 5.0e-4);
  return 0;
}
