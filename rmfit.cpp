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
  int nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);

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
                << "blocks = ...\n"
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
  
  Log::open(params.output.c_str());
  
  if (Log::ionode)
  {
    log() << "# RMFIT MPI v4.0\n"
          << "Parameters provided:\n"
          << "  N      = " << params.N << '\n'
          << "  nu     = " << params.nu << '\n'
          << "  sigma  = " << params.center.coord[0];
    if (params.active[0])
      log() << '\t' << params.scale.coord[0];
    log() << "\n  m      = " << params.center.coord[1];
    if (params.active[1])
      log() << '\t' << params.scale.coord[1];
    log() << "\n  a6     = " << params.center.coord[2];
    if (params.active[2])
      log() << '\t' << params.scale.coord[2];
    log() << "\n  a7     = " << params.center.coord[3];
    if (params.active[3])
      log() << '\t' << params.scale.coord[3];
    log() << "\n  a8     = " << params.center.coord[4];
    if (params.active[4])
      log() << '\t' << params.scale.coord[4];
    log() << "\n  blocks = " << params.blocks << '\n'
          << "  prec   = " << params.prec << '\n'
          << "  data   = " << params.data << '\n'
          << "  output = " << params.output << std::endl;
    if (params.dim == 1)
      log() << "There are no degrees of freedom in the fit,\n"
            << "so this run will just determine the results\n"
            << "for the parameters given above." << std::endl;
    log() << "Run on " << nodes << " nodes." << std::endl;
  }
  
  Data data(params.data.c_str());
  Minim minim(data, params);

  Simplex const &result = minim.reduce();
  
  Log::shut();
  
  MPI_Finalize();
  return 0;
}
