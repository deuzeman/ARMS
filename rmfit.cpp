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
      std::cerr << "# RMFIT MPI v3.1\n"
		<< "Need an input file!\n" 
		<< "This file should contain the following lines (... indicate input values):\n\n"
		<< "data = ...\n"
		<< "N = ...\n"
		<< "nu = ...\n"
		<< "m = ... ...\n"
		<< "a6 = ... ...\n"
		<< "a7 = ... ...\n"
		<< "a8 = ... ...\n"
		<< "sigma = ... ...\n"
		<< "iter = ... ...\n"
                << "tol = ...\n"
                << "kolmogorov = ...\n"
		<< "bootSeed = ..\n\n"
		<< "The name of this input file should be the only argument.\n"
		<< "Now exiting." << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  FitParams params(argv[1], true); // Reads in parallel
  Data data(params.data.c_str(), params.bootSeed, true); // Reads in parallel
   
  Minim *minim = new Minim(data, params.N, params.nu, params.kol); // Set up the minimization structure

  if (myRank != 0)
  {
    minim->listen();
  }
  else
  {
    int mpisize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    // If we arrive here, we're the master node.
    minim->log() << "**********************************\n"
                << "RMFIT MPI v3.1 -- FITTING REPORT\n\n"
                << "data = " << params.data << '\n'
                << "N = " << params.N << '\n'
                << "nu = " << params.nu << '\n'
                << "m = " << params.m[0] <<  ' ' << params.m[1] << '\n'
                << "a6 = " << params.a6[0] << ' ' << params.a6[1] << '\n'
                << "a7 = " << params.a7[0] << ' ' << params.a7[1] << '\n'
                << "a8 = " << params.a8[0] << ' ' << params.a8[1] << '\n'
                << "sigma = " << params.sigma[0] << ' ' << params.sigma[1] << '\n'
                << "iter = " << params.iter[0] << ' ' << params.iter[1] << '\n'
                << "tol = " << params.tol << '\n'
                << "kolmogorov = " << params.kol << '\n'
                << "bootSeed = " << params.bootSeed << '\n'
                << "Running with " << mpisize << " MPI processes\n"
                << "**********************************\n" << std::endl;
    
    // We want to jiggle the starting values a bit, to avoid sklerotic behaviour
    // The quality of this doesn't have to be brilliant, so let's just use C's crappy random function
    srand(time(0));
    Eigen::ArrayXd start(5);
    start << params.m[0] + (params.m[1] - params.m[0]) * (static_cast< double >(rand()) / RAND_MAX),
             params.a6[0] + (params.a6[1] - params.a6[0]) * (static_cast< double >(rand()) / RAND_MAX),
             params.a7[0] + (params.a7[1] - params.a7[0]) * (static_cast< double >(rand()) / RAND_MAX),
             params.a8[0] + (params.a8[1] - params.a8[0]) * (static_cast< double >(rand()) / RAND_MAX),
             params.sigma[0] + (params.sigma[1] - params.sigma[0]) * (static_cast< double >(rand()) / RAND_MAX);

    Eigen::ArrayXXd bounds(2, 5);
    bounds << params.m[0], params.a6[0], params.a7[0], params.a8[0], params.sigma[0],
              params.m[1], params.a6[1], params.a7[1], params.a8[1], params.sigma[1];
    minim->powell(start, bounds, params.iter[0], params.iter[1], 50, params.tol);
  }
  delete minim;
  MPI_Finalize();
  return 0;
}
