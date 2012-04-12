#include <fstream>
#include <mpi.h>
#include <sstream>

#include <Params.h>
#include <Log.h>

// This is as good a point as any to set up the mapping in point!
//   coord[0] = sigma
//   coord[1] = m
//   coord[2] = a6
//   coord[3] = a7
//   coord[4] = a8

static void cleanLine(std::string &line, char const *buff)
{
  char clBuff[256];
  char *cpb = clBuff;
  
  bool before = true;
  bool clear = true;
  for (char const *pb = buff; *pb != '\0'; ++pb)
  {
    if (clear && ((*pb == '\t' ) || (*pb == '\n' ) || (*pb == ' ' )))
      continue;

    if (!before)
      clear = false;
    
    if (*pb == '=')
      before = false;
    
    *cpb = *pb;
    ++cpb;
  }
  *cpb = '\0';
  line.assign(clBuff);
}

Params::Params(char const *filename)
{
  int nodes;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::ifstream input;
  char line[256];
  std::string sline;
  
  if (rank == 0)
  {
    input.open(filename, std::ifstream::in);
  
    while(input.getline(line, 256))
    {
       cleanLine(sline, line);

      if (sline.find("N=") != sline.npos)
      {
        N = XLat(sline.c_str() + 2);
        continue;
      }

      if (sline.find("nu=") != sline.npos)
      {
        nu = XLat(sline.c_str() + 3);
        continue;
      }

      if (sline.find("sigma=") != sline.npos)
      {
        XLat sigma(sline.c_str() + 6);
        center.coord[0] = sigma;
        scale.coord[0] = sigma;
        continue;
      }

      if (sline.find("m=") != sline.npos)
      {
        XLat mass(sline.c_str() + 2);
        center.coord[1] = mass;
        scale.coord[1] = mass;
        continue;
      }
      
      if (sline.find("a6=") != sline.npos)
      {
        XLat a6(sline.c_str() + 3);
        center.coord[2] = a6;
        scale.coord[2] = a6;
        continue;
      }
      
      if (sline.find("a7=") != sline.npos)
      {
        XLat a7(sline.c_str() +3);
        center.coord[3] = a7;
        scale.coord[3] = a7;
        continue;
      }

      if (sline.find("a8=") != sline.npos)
      {
        XLat a8(sline.c_str() + 3);
        center.coord[4] = a8;
        scale.coord[4] = a8;
        continue;
      }
      
      if (sline.find("data=") != sline.npos)
      {
        data.assign(sline.c_str() + 5);
        continue;
      }

      if (sline.find("output=") != sline.npos)
      {
        output.assign(sline.c_str() + 7);
        continue;
      }
      
      if (sline.find("prec=") != sline.npos)
      {
        prec = XLat(sline.c_str() + 5);
        continue;
      }      
     
      if (sline.find("eigMin=") != sline.npos)
      {
        eigMin = XLat(sline.c_str() + 7);
        continue;
      }
      
      if (sline.find("eigMax=") != sline.npos)
      {
        eigMax = XLat(sline.c_str() + 7);
        continue;
      }
      
      if (sline.find("iter=") != sline.npos)
      {
        iter = XLat(sline.c_str() + 5);
        continue;
      }  
    }
  }

  // Now work through the whole parameters struct and broadcast everything
  if (rank == 0)
  {
    std::copy(data.c_str(), data.c_str() + data.length(), line);
    line[data.length()] = '\0';
  }
  MPI_Bcast(line, data.length() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (rank != 0)
    data.assign(line);

  if (rank == 0)
  {
    std::copy(output.c_str(), output.c_str() + output.length(), line);
    line[output.length()] = '\0';
  }
  MPI_Bcast(line, output.length() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (rank != 0)
    output.assign(line);
  
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nu, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(center.coord, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(scale.coord, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&prec, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eigMin, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eigMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
}