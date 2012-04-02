#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <Params.h>

void Params::parseInput(char const *filename)
{
  std::ifstream input(filename);
  char line[256];

  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      m = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a6 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a7 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a8 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("sigma = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      sigma = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nEig_min = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEig_min = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nEig_max = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEig_max = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nDet = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nDet = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      iter = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("output = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      output = (sline.c_str() + idx + 2);
      continue;
    }
  }
}

void FitParams::parseInput(char const *filename)
{
  std::ifstream input(filename);
  char line[256];

  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("data = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      std::stringstream sst(sline.c_str() + idx + 2);
      sst >> data;
      continue;
    }

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      m[0] = istr;
      m[1] = istr;
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a6[0] = istr;
      a6[1] = istr;
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a7[0] = istr;
      a7[1] = istr;
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a8[0] = istr;
      a8[1] = istr;
      continue;
    }

    if (sline.find("sigma = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      sigma[0] = istr;
      sigma[1] = istr;
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      iter[0] = istr;
      iter[1] = istr;
      continue;
    }

    if (sline.find("tol = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      tol = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("kolmogorov = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      kol = XLat(sline.c_str() + idx + 2);
      continue;
    }
    
    if (sline.find("bootSeed = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      bootSeed = XLat(sline.c_str() + idx + 2);
      continue;
    }
  }
}

void FitParams::parseInputParallel(char const *filename)
{
  int mpirank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  // Broken LUSTRE forces us to use POSIX, then broadcast
  unsigned long int filesize = 0;
  if (mpirank == 0)
  {
    struct stat filestatus;
    stat(filename, &filestatus);
    filesize = filestatus.st_size;
  }
  MPI_Bcast(&filesize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  char *buffer = new char[filesize + 1]; // To avoid weirdness, assign a little extra memory

  if (mpirank == 0)
  {
    std::ifstream instream(filename);
    instream.read(buffer, filesize);
    instream.close();
  }

  MPI_Bcast(buffer, filesize, MPI_CHAR, 0, MPI_COMM_WORLD);
  buffer[filesize] = '\0'; // For security reasons, don't know if it's really needed, can't hurt.

  std::istringstream input(buffer);
  char line[256];
  // Should be business as usual from here on!
  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("data = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      std::stringstream sst(sline.c_str() + idx + 2);
      sst >> data;
      continue;
    }

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      m[0] = istr;
      m[1] = istr;
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a6[0] = istr;
      a6[1] = istr;
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a7[0] = istr;
      a7[1] = istr;
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a8[0] = istr;
      a8[1] = istr;
      continue;
    }

    if (sline.find("sigma = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      sigma[0] = istr;
      sigma[1] = istr;
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      iter[0] = istr;
      iter[1] = istr;
      continue;
    }

    if (sline.find("tol = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      tol = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("kolmogorov = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      kol = XLat(sline.c_str() + idx + 2);
      continue;
    }
    
    if (sline.find("bootSeed = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      bootSeed = XLat(sline.c_str() + idx + 2);
      continue;
    }
  }
  delete[] buffer;
}
