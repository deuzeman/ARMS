#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <RanMat.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

int main(int argc, char **argv)
{
  if (argc == 1)
  {
    std::cerr << "# ARMS v2.0\n"
              << "Need an input file!\n" 
              << "This file should contain the following lines:\n\n"
              << "N = ...\n"
              << "nu = ...\n"
              << "m = ...\n"
              << "a6 = ...\n"
	      << "a7 = ...\n"
	      << "a8 = ...\n"
              << "iter = ...\n"
              << "output = ...\n\n"
              << "The name of this input file should be the first argument.\n"
              << "Now exiting." << std::endl;
    return 1;
  }

  params const *par = parseInput(argv[1]);
  RanMat generator(par, time(0));

  std::ofstream rstream(par->output.c_str(), std::ifstream::trunc);
  rstream << "# ARMS v2.0 output file\n"
          << "# Scaling massless operator with inverse sqrt(N)\n"
          << "# Extended version, with three Wilson operators in the action.\n"
          << "# Run with the following parameters\n"
          << "# N:    " << par->N    << '\n'
          << "# nu:   " << par->nu   << '\n'
          << "# m:    " << par->m    << '\n'
          << "# a6:   " << par->a6    << '\n'
          << "# a7:   " << par->a7    << '\n'
          << "# a8:   " << par->a8    << '\n'
          << "# iter: " << par->iter << '\n'
          << '#'                     << std::endl;

  for (int x = -6; x < 6; ++x)
  {
    int colNum = (x < 0) ? x : x + 1;
    if (x < 0)
      rstream << "\"EV.m" << -colNum << "\" ";
    else
      rstream << "\"EV.p" << colNum << "\" ";
  }
  rstream << std::endl;

  generator.calculate(par->iter);

  for (int k = 0; k < par->iter; ++k)
  {
    rstream << "\"" << (k + 1) << "\" ";

    for (int x = 0; x < 12; ++x)
      rstream << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(8) << generator.result(k, x);
    rstream << std::endl;
  }
  rstream.close();
  delete par;

  return 0;
}
