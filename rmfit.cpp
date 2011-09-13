#include <iostream>
#include <cstdlib>

#include <Data.h>
#include <Minim.h>
#include <Params.h>

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cerr << "# RMFIT v2.3\n"
              << "Need an input file!\n" 
              << "This file should contain the following lines (... indicate input values):\n\n"
              << "data = ...\n"
              << "N = ...\n"
              << "nu = ...\n"
              << "m = ... ... ...\n"
              << "a6 = ... ... ...\n"
              << "a7 = ... ... ...\n"
              << "a8 = ... ... ...\n"
              << "scale = ... ... ...\n"
              << "iter = ... ...\n"
              << "tol = ...\n\n"
              << "The name of this input file should be the only argument.\n"
              << "Now exiting." << std::endl;
    return 1;
  }
  FitParams params(argv[1]);
  Data data(params.data.c_str());
  Minim minim(data, params.N, params.nu, data.minEv(), data.maxEv(), time(0), true);

  minim.log() << "RMFIT v2.3 -- FITTING REPORT\n"
              << "Called with the following options:\n"
              << "data:     " << params.data << '\n'
              << "N:        " << params.N << '\n'
              << "nu:       " << params.nu << '\n'
              << "m:        " << params.m[0] << ' ' << params.m[1] << ' ' << params.m[2] << '\n'
              << "a6:       " << params.a6[0] << ' ' << params.a6[1] << ' ' << params.a6[2] << '\n'
              << "a7:       " << params.a7[0] << ' ' << params.a7[1] << ' ' << params.a7[2] << '\n'
              << "a8:       " << params.a8[0] << ' ' << params.a8[1] << ' ' << params.a8[2] << '\n'
              << "scale:    " << params.scale[0] << ' ' << params.scale[1] << ' ' << params.scale[2] << '\n'
              << "iter:     " << params.iter[0] << ' ' << params.iter[1] << '\n'
              << "tol:      " << params.tol << '\n' << std::endl;

  Eigen::ArrayXd start(5);
  start << params.m[1], params.a6[1], params.a7[1], params.a8[1], params.scale[1];
  Eigen::ArrayXXd bounds(2, 5);
  bounds << params.m[0], params.a6[0], params.a7[0], params.a8[0], params.scale[0],
            params.m[2], params.a6[2], params.a7[2], params.a8[2], params.scale[2];
  minim.powell(start, bounds, params.iter[0], params.iter[1], 50, params.tol);

  return 0;
}
