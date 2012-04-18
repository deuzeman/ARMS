#pragma once

#include <mpi.h>

#include <string>
#include <sstream>
#include <Point.h>

class XLat
{
  std::istringstream d_data;

  public:
    XLat(char const *input)
      : d_data(input)
    {};

    template< typename Type >
    operator Type()
    {
      Type result;
      d_data >> result;
      return result;
    }
    
    bool good() const
    {
      return d_data.good();
    }
};

struct Params
{
  std::string data;
  std::string output;

  size_t N;
  size_t nu;

  Point center;
  Point scale;
  
  double prec;
  size_t blocks;
  
  int eigMin;
  int eigMax;
  
  bool active[5];
  unsigned int dim;
  unsigned int iter;

  Params(char const *filename);
};
