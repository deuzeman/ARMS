#include <iostream>
#include <cstdlib>

#include <mpi.h>

#include <Data.h>
#include <Minim.h>
#include <Params.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if (argc != 2)
  {
    if (myRank == 0)
    {
      std::cerr << "# RMFIT MPI v4.0\n"
		<< "Need an input file!\n" 
		<< "This file should contain the following lines (... indicate input values):\n\n"
		<< "N      = ...\n"
		<< "nu     = ...\n"
                << "sigma  = ... ...\n"
                << "m      = ... ...\n"
                << "a6     = ... ...\n"
                << "a7     = ... ...\n"
                << "a8     = ... ...\n"
                << "prec   = ...\n"
                << "data   = ...\n"
                << "output = ...\n\n"
                << "Double arguments for the parameters give the best guess (first) and scale (second)\n"
		<< "The name of this input file should be the only argument.\n"
		<< "Now exiting." << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  Params params(argv[1]);
  
  if (myRank == 0)
  {
    Log::open(params.output.c_str());
    log() << "# RMFIT MPI v4.0\n"
          << "Parameters provided:\n"
          << "  N      = " << params.N << '\n'
          << "  nu     = " << params.nu << '\n'
          << "  sigma  = " << params.center.coord[0] << '\t' << params.scale.coord[0] << '\n'
          << "  m      = " << params.center.coord[1] << '\t' << params.scale.coord[1] << '\n'
          << "  a6     = " << params.center.coord[2] << '\t' << params.scale.coord[2] << '\n'
          << "  a7     = " << params.center.coord[3] << '\t' << params.scale.coord[3] << '\n'
          << "  a8     = " << params.center.coord[4] << '\t' << params.scale.coord[4] << '\n'
          << "  prec   = " << params.prec << '\n'
          << "  data   = " << params.data << '\n'
          << "  output = " << params.output << '\n' << std::endl;
  }
  
  Data data(params.data.c_str());
  Minim minim(data, params);

  Simplex const &result = minim.reduce();
  
  if (myRank == 0)
    Log::shut();
  
  MPI_Finalize();
  return 0;
}
