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
    std::cerr << "# RMZERO v1.1\n"
              << "Need an input file!\n" 
              << "This file should contain the following lines:\n\n"
              << "N = ...\n"
              << "nu = ...\n"
              << "m = ...\n"
              << "a6 = ...\n"
              << "a7 = ...\n"
              << "a8 = ...\n"
              << "sigma = ...\n"
              << "iter = ...\n"
              << "output = ...\n\n"
              << "The name of this input file should be the only argument.\n"
              << "Now exiting." << std::endl;
    return 1;
  }

  Params par(argv[1]);
  RanMat generator(par.N, par.nu, -par.N, +par.N, 0, time(0), true);

  std::ofstream rstream(par.output.c_str(), std::ofstream::trunc);
  rstream << "# RMZERO v1.1 output file\n"
          << "# Run with the following parameters\n"
          << "# N:        " << par.N        << '\n'
          << "# nu:       " << par.nu       << '\n'
          << "# m:        " << par.m        << '\n'
          << "# a6:       " << par.a6       << '\n'
          << "# a7:       " << par.a7       << '\n'
          << "# a8:       " << par.a8       << '\n'
          << "# sigma:    " << par.sigma    << '\n'
          << "# iter:     " << par.iter     << '\n'
          << '#'                            << std::endl;


  Eigen::ArrayXd rmtPars(5);
  rmtPars << par.m, par.a6, par.a7, par.a8, par.sigma;

  for (int k = 0; k < par.iter; ++k)
  {
    generator.calcreal(rmtPars);

    Eigen::Matrix< std::complex< double >, Eigen::Dynamic, 1 > const &ev = generator.eigenvalues_cmplx();
    Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > const & vecs = generator.eigenvectors_cmplx();
    for (int x = 0; x < ev.rows(); ++x)
    {
      if (std::abs(std::imag(ev[x])) < 1E-10)
      {
        rstream << std::setw(10) << std::setprecision(8) << std::real(ev[x])
                << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(8) << std::real((vecs.col(x).adjoint() * generator.gamma_5() * vecs.col(x)).eval()[0]) << std::endl;
      }
    }
  }
  rstream.close();

  std::cerr << "Done." << std::endl;

  return 0;
}
