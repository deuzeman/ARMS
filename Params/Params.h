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
  double sigma;
  int    nEig_min;
  int    nEig_max;
  size_t nDet;
  size_t iter;
  std::string output;

  Params(char const *filename);
  void parseInput(char const *filename);
};

inline Params::Params(char const *filename)
{
  parseInput(filename);
}

struct FitParams
{
  std::string data;

  size_t N;
  size_t nu;

  double m[2];
  double a6[2];
  double a7[2];
  double a8[2];
  double sigma[2];

  size_t iter[2];
  double tol;
  size_t kol;
  int bootSeed;

  FitParams(char const *filename, bool const parallel = false);
  void parseInput(char const *filename);
  void parseInputParallel(char const *filename);
};

inline FitParams::FitParams(char const *filename, bool const parallel)
  : kol(0)
{
  if (parallel)
    parseInputParallel(filename);
  else
   parseInput(filename);
}
