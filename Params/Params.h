#pragma once

#include <string>
#include <sstream>

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
};

struct Params
{
  size_t N;
  size_t nu;
  double m;
  double a6;
  double a7;
  double a8;
  int    nEigs;
  size_t nDet;
  size_t iter;
  std::string output;
};

Params *parseInput(char const *filename);
