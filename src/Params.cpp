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
  
  // Add code to trim the string at the end, or istringstream will not set the EOF flag correctly
  char const *pbl = buff;
  while (*pbl != '\0')
    ++pbl;
  while ((*pbl == '\0') || (*pbl == '\t' ) || (*pbl == '\n' ) || (*pbl == ' '))
    --pbl;
  ++pbl;
  
  bool before = true;
  bool clear = true;
  for (char const *pb = buff; pb != pbl; ++pb)
  {
    bool white = ((*pb == '\t' ) || (*pb == '\n' ) || (*pb == ' ' ));
    
    if (clear && white)
      continue;

    if (!before)
      clear = white; // Remove if subsequent whitespace!
      
    if (*pb == '=')
    {
      before = false;
      clear = true;
    }
       
    *cpb = *pb;
    ++cpb;
  }
  *cpb = '\0';
  line.assign(clBuff);
}

Params::Params(char const *filename)
  : dim(1)
{
  int nodes;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::ifstream input;
  char line[256];
  std::string sline;
  std::fill_n(active, 5, false);
  std::fill_n(scale.coord, 5, 0.0);
  
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
        if (sigma.good())
        {
          ++dim;
          active[0] = true;
          scale.coord[0] = sigma;
        }
        continue;
      }

      if (sline.find("m=") != sline.npos)
      {
        XLat mass(sline.c_str() + 2);
        center.coord[1] = mass;
        if (mass.good())
        {
          ++dim;
          active[1] = true;
          scale.coord[1] = mass;
        }
        continue;
      }
      
      if (sline.find("a6=") != sline.npos)
      {
        XLat a6(sline.c_str() + 3);
        center.coord[2] = a6;
        if (a6.good())
        {
          ++dim;
          active[2] = true;
          scale.coord[2] = a6;
        }
        continue;
      }
      
      if (sline.find("a7=") != sline.npos)
      {
        XLat a7(sline.c_str() +3);
        center.coord[3] = a7;
        if (a7.good())
        {
          ++dim;
          active[3] = true;
          scale.coord[3] = a7;
        }
        continue;
      }

      if (sline.find("a8=") != sline.npos)
      {
        XLat a8(sline.c_str() + 3);
        center.coord[4] = a8;
        if (a8.good())
        {
          ++dim;
          active[4] = true;
          scale.coord[4] = a8;
        }
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
      
      if (sline.find("blocks=") != sline.npos)
      {
        blocks = XLat(sline.c_str() + 7);
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

  int toCommunicate = 0;
  if (rank == 0)
    toCommunicate = data.length();
  MPI_Bcast(&toCommunicate, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Now work through the whole parameters struct and broadcast everything
  if (toCommunicate > 0)
  {
    if (rank == 0)
    {
      std::copy(data.c_str(), data.c_str() + toCommunicate, line);
      line[toCommunicate] = '\0';
    }
    MPI_Bcast(line, toCommunicate + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0)
      data.assign(line);
  }
  
  if (rank == 0)
    toCommunicate = output.length();
  MPI_Bcast(&toCommunicate, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (toCommunicate > 0)
  {
    if (rank == 0)
    {
      std::copy(output.c_str(), output.c_str() + toCommunicate, line);
      line[toCommunicate] = '\0';
    }
    MPI_Bcast(line, toCommunicate + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0)
      output.assign(line);
  }
  
  int temp = static_cast< int >(N);
  MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  N = static_cast< unsigned int >(temp);

  temp = static_cast< int >(dim);
  MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dim = static_cast< unsigned int >(temp);
  
  temp = static_cast< int >(nu);
  MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  nu = static_cast< unsigned int >(temp);
  
  temp = static_cast< int >(iter);
  MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  iter = static_cast< unsigned int >(temp);
  
  temp = static_cast< int >(blocks);
  MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  blocks = static_cast< unsigned int >(temp);
  
  MPI_Bcast(center.coord, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(scale.coord, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&prec, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eigMin, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eigMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
}