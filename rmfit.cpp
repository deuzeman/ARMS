#include <Data.h>

int main(int argc, char **argv)
{
  Data data("test.dat");
  std::cout << data.average(1000).row(1) / data.average(1000).row(0) << std::endl;
  std::cout << data.ratios(1000).row(1) / data.ratios(1000).row(0) << std::endl;
  std::cout << data.cumulant().leftCols(3) << std::endl;
  return 0;
}
