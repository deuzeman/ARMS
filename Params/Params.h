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
  double scale;
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

  double m[3];
  double a6[3];
  double a7[3];
  double a8[3];
  double scale[3];

  size_t iter[2];
  double tol;

  FitParams(char const *filename);
  void parseInput(char const *filename);
};

inline FitParams::FitParams(char const *filename)
{
  parseInput(filename);
}