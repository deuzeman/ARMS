#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <mpi.h>

#include <Data.h>
#include <Log.h>
#include <Minim.h>
#include <Params.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  int nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);
  
  if (argc != 2)
  {
    if (myRank == 0)
    {
      std::cerr << "# RMS MPI v4.0\n"
      << "Need an input file!\n" 
      << "This file should contain the following lines (... indicate input values):\n\n"
      << "N      = ...\n"
      << "nu     = ...\n"
      << "sigma  = ...\n"
      << "m      = ...\n"
      << "a6     = ...\n"
      << "a7     = ...\n"
      << "a8     = ...\n"
      << "eigMin = ...\n"
      << "eigMax = ...\n"
      << "iter   = ...\n"
      << "output = ...\n\n"
      << "Double arguments for the parameters give the best guess (first) and scale (second)\n"
      << "The name of this input file should be the only argument.\n"
      << "Now exiting." << std::endl;
    }
    MPI_Finalize();
    return 1;
  }
  
  Params params(argv[1]);
  
  Log::open(params.output.c_str());
  
  if (myRank == 0)
  {
    log() << "# RMS MPI v4.0\n"
          << "# Parameters provided:\n"
          << "#   N      = " << params.N << '\n'
          << "#   nu     = " << params.nu << '\n'
          << "#   sigma  = " << params.center.coord[0] << '\n'
          << "#   m      = " << params.center.coord[1] << '\n'
          << "#   a6     = " << params.center.coord[2] << '\n'
          << "#   a7     = " << params.center.coord[3] << '\n'
          << "#   a8     = " << params.center.coord[4] << '\n'
          << "#   eigMin = " << params.eigMin << '\n'
          << "#   eigMax = " << params.eigMax << '\n'
          << "#   iter   = " << params.iter << '\n'
          << "#   output = " << params.output << '\n' 
          << "# Run on " << nodes << " nodes." << '\n'
          << "#"             << std::endl;
  }
  
  RanMat generator(params.N, params.nu, params.eigMin, params.eigMax);
  
  generator.calculate(params.center, params.iter);
  double const *result = generator.result();
  
  if (Log::ionode)
  {
    log().setf(std::ios::right, std::ios::adjustfield);
    for (int ctr = params.eigMin; ctr <= params.eigMax; ++ctr)
    {
      if (!ctr)
        continue;
      log() << std::setw(14) << (ctr < 0 ? 'M' : 'P')  << std::abs(ctr);
    }
    log() << std::endl;
    Log::shut();
  }
  
  for (size_t idx = 0; idx < generator.nodes(); ++idx)
  {
    Log::open(params.output.c_str(), idx);
    if (Log::ionode)
    {
      log().precision(8);
      log().setf(std::ios::fixed, std::ios::floatfield);
      log().setf(std::ios::right, std::ios::adjustfield);
      for (size_t it = 0; it < generator.numSamples(); ++it)
      {
        for (size_t eig = 0; eig < generator.numEigs(); ++eig)
          log() << std::setw(15) << result[eig * generator.numSamples() + it];
        log() << std::endl;
      }
    }
    Log::shut();
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Finalize();
  return 0;
}

