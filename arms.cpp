#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

#include <Params.h>
#include <RanMat.h>

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cerr << "# ARMS v2.5\n"
              << "Need an input file!\n" 
              << "This file should contain the following lines:\n\n"
              << "N = ...\n"
              << "nu = ...\n"
              << "m = ...\n"
              << "a6 = ...\n"
              << "a7 = ...\n"
              << "a8 = ...\n"
              << "sigma = ...\n"
              << "nEig_min = ...\n"
              << "nEig_max = ...\n"
              << "nDet = ...\n"
              << "iter = ...\n"
              << "output = ...\n\n"
              << "The name of this input file should be the only argument.\n"
              << "Now exiting." << std::endl;
    return 1;
  }

  Params par(argv[1]);
  RanMat generator(par.N, par.nu, par.nEig_min, par.nEig_max, par.nDet, time(0), true);

  std::ofstream rstream(par.output.c_str(), std::ofstream::trunc);
  rstream << "# ARMS v2.5 output file\n"
          << "# Scaling massless operator with inverse sqrt(N)\n"
          << "# Extended version, with three Wilson operators in the action.\n"
          << "# Run with the following parameters\n"
          << "# N:        " << par.N        << '\n'
          << "# nu:       " << par.nu       << '\n'
          << "# m:        " << par.m        << '\n'
          << "# a6:       " << par.a6       << '\n'
          << "# a7:       " << par.a7       << '\n'
          << "# a8:       " << par.a8       << '\n'
          << "# sigma:    " << par.sigma    << '\n'
          << "# nEig_min: " << par.nEig_min << '\n'
          << "# nEig_max: " << par.nEig_max << '\n'
          << "# nDet:     " << par.nDet     << '\n'
          << "# iter:     " << par.iter     << '\n'
          << '#'                            << std::endl;

  rstream << std::setw(10) << ' ';
  for (int x = par.nEig_min; x < par.nEig_max; ++x)
  {
    int colNum = (x < 0) ? x : x + 1;
    std::ostringstream colLab;
    if (x < 0)
      colLab << "\"EV.m" << -colNum << '\"';
    else
      colLab << "\"EV.p" << colNum << '\"';
   rstream << std::setw(15) << colLab.str();
  }
  if (par.nDet > 0)
    rstream << std::setw(15) << "\"Det\"";
  rstream << std::endl;

  Eigen::ArrayXd rmtPars(5);
  rmtPars << par.m, par.a6, par.a7, par.a8, par.sigma;
  generator.calculate(rmtPars, par.iter, false);

  for (int k = 0; k < par.iter; ++k)
  {
    std::ostringstream lineLab;
    lineLab << '\"' << (k + 1) << '\"';
    rstream << std::setw(10) << lineLab.str();

    for (int x = 0; x < par.nEig_max - par.nEig_min; ++x)
      rstream << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(8) << generator.result(k, x);

    if (par.nDet > 0)
      rstream << std::setw(15) << std::setiosflags(std::ios::scientific) << generator.det(k);

    rstream << std::endl;
  }
  rstream.close();

  std::cerr << "Done. Averages: " << std::endl;
  std::cerr << generator.average() << std::endl;

  return 0;
}
