#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

#include <RanMat.h>
#include <Params.h>

int main(int argc, char **argv)
{
  if (argc == 1)
  {
    std::cerr << "# ARMS v2.1\n"
              << "Need an input file!\n" 
              << "This file should contain the following lines:\n\n"
              << "N = ...\n"
              << "nu = ...\n"
              << "m = ...\n"
              << "a6 = ...\n"
              << "a7 = ...\n"
              << "a8 = ...\n"
              << "nEigs = ...\n"
              << "nDet = ...\n"
              << "iter = ...\n"
              << "output = ...\n\n"
              << "The name of this input file should be the first argument.\n"
              << "Now exiting." << std::endl;
    return 1;
  }

  Params *par = parseInput(argv[1]);
  RanMat generator(par->N, par->nu, par->nEigs, par->nDet, time(0));

  std::ofstream rstream(par->output.c_str(), std::ifstream::trunc);
  rstream << "# ARMS v2.0 output file\n"
          << "# Scaling massless operator with inverse sqrt(N)\n"
          << "# Extended version, with three Wilson operators in the action.\n"
          << "# Run with the following parameters\n"
          << "# N:     " << par->N     << '\n'
          << "# nu:    " << par->nu    << '\n'
          << "# m:     " << par->m     << '\n'
          << "# a6:    " << par->a6    << '\n'
          << "# a7:    " << par->a7    << '\n'
          << "# a8:    " << par->a8    << '\n'
          << "# nEigs: " << par->nEigs << '\n'
          << "# nDet:  " << par->nDet  << '\n'
          << "# iter:  " << par->iter  << '\n'
          << '#'                       << std::endl;

  rstream << "                  ";
  for (int x = -(par->nEigs / 2); x < (par->nEigs / 2 + par->nEigs % 2); ++x) {
    int colNum = (x < 0) ? x : x + 1;
    if (x < 0)
      rstream << "\"EV.m" << -colNum << "\"        ";
    else
      rstream << "\"EV.p" <<  colNum << "\"        ";
  }
  if (par->nDet > 0)
    rstream << "\"Det\"";
  rstream << std::endl;

  Eigen::ArrayXd rmtPars(4);
  rmtPars << par->m, par->a6, par->a7, par->a8;
  generator.calculate(rmtPars, par->iter, false);

  for (int k = 0; k < par->iter; ++k)
  {
    std::ostringstream lineLab;
    lineLab << '\"' << (k + 1) << '\"';
    rstream << std::setw(10) << lineLab.str();

    for (int x = 0; x < par->nEigs; ++x)
      rstream << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(8) << generator.result(k, x);

    if (par->nDet > 0)
      rstream << std::setw(15) << std::setiosflags(std::ios::scientific) << generator.det(k);

    rstream << std::endl;
  }
  rstream.close();
  delete par;

  return 0;
}
