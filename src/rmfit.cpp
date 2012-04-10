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
		<< "data = ...\n"
		<< "N = ...\n"
		<< "nu = ...\n"
		<< "m = ...\n"
		<< "a6 = ...\n"
		<< "a7 = ...\n"
		<< "a8 = ...\n"
		<< "sigma = ...\n"
		<< "The name of this input file should be the only argument.\n"
		<< "Now exiting." << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  FitParams params(argv[1], true); // Reads in parallel
  Data data(params.data.c_str(), params.bootSeed, true); // Reads in parallel
   
  Minim *minim = new Minim(&data, &params);
  Simplex const &result = minim->reduce();
  
  if (myRank == 0)
  {
    std::cout << "Obtained the following result:\n"
              << "Values for D: " << result.values[0].value << ' ' << result.values[1].value << ' ' 
                                  << result.values[2].value << ' ' << result.values[3].value << ' '
                                  << result.values[4].value << ' ' << result.values[5].value << '\n'
              << "Values for sigma: " << result.points[0].sigma << ' ' << result.points[1].sigma << ' ' 
                                  << result.points[2].sigma << ' ' << result.points[3].sigma << ' '
                                  << result.points[4].sigma << ' ' << result.points[5].sigma << '\n'
              << "Values for m: " << result.points[0].m << ' ' << result.points[1].m << ' ' 
                                  << result.points[2].m << ' ' << result.points[3].m << ' '
                                  << result.points[4].m << ' ' << result.points[5].m << '\n'
              << "Values for a6: " << result.points[0].a6 << ' ' << result.points[1].a6 << ' ' 
                                  << result.points[2].a6 << ' ' << result.points[3].a6 << ' '
                                  << result.points[4].a6 << ' ' << result.points[5].a6 << '\n'
              << "Values for a7: " << result.points[0].a7 << ' ' << result.points[1].a7 << ' ' 
                                  << result.points[2].a7 << ' ' << result.points[3].a7 << ' '
                                  << result.points[4].a7 << ' ' << result.points[5].a7 << '\n'
              << "Values for a8: " << result.points[0].a8 << ' ' << result.points[1].a8 << ' ' 
                                  << result.points[2].a8 << ' ' << result.points[3].a8 << ' '
  }
  delete minim;
  MPI_Finalize();
  return 0;
}
