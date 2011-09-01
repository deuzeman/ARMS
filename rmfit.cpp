#include <iostream>
#include <Data.h>
#include <Minim.h>

int main(int argc, char **argv)
{
  Data data("test.dat");
  Minim minim(data, 20, 1, data.minEv(), data.maxEv(), 234231);
  Eigen::ArrayXd start(2);
  start << 0.05, 0.02;
  Eigen::ArrayXXd bounds(2,2);
  bounds << 0.0, 0.0, 0.1, 0.05; 
  minim.powell(start, bounds, 1000, 50, 5e-3);
  return 0;
}
