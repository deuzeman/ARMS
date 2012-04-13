#include <Log.h>

#include <mpi.h>
#include <cstdlib>

Log *Log::s_instance = 0;
bool Log::ionode = false;

void Log::open(char const *filename, int rank)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ionode = (myrank == rank);
  if (ionode)
  {
    if (Log::s_instance)
    {
      std::cerr << "[ERROR] Double opening of log!" << std::endl;
      exit(1);
    }
    Log::s_instance = new Log(filename);
  }
}

void Log::shut()
{
  if (s_instance)
    s_instance->close();
  delete Log::s_instance;
  Log::s_instance = 0;
  Log::ionode = false;
}

std::ofstream &Log::put()
{
  if (!Log::s_instance)
  {
    std::cerr << "[ERROR] Attempted to use log without opening!" << std::endl;
    exit(1);
  }
  return *s_instance;
}