#pragma once

#include <cstdlib>
#include <fstream>
#include <sstream>

#include <Eigen/Dense>
#include <Random/sfmt.h>


class Data
{
  int d_minEv;
  int d_maxEv;
  size_t d_nSamp;

  Eigen::ArrayXXd d_data;
  double d_normalization;

  Eigen::ArrayXd  d_average;
  Eigen::ArrayXd  d_error;
  Eigen::ArrayXXd d_cumulant;

  public:
    Data(char const *filename);

    double *average() const;
    double *error() const;
    Eigen::ArrayXXd const &cumulant() const;

    int minEv() const;
    int maxEv() const;
    size_t numSamples() const;
    size_t numCols() const;
    int bootSeed() const;

    double normalization() const;

    double *flat() const;
    double *flat(size_t const &col) const;
    double *flatPerColumn() const;
};

inline int Data::minEv() const
{
  return d_minEv;
}

inline int Data::maxEv() const
{
  return d_maxEv;
}

inline size_t Data::numSamples() const
{
  return d_nSamp;
}

inline size_t Data::numCols() const
{
  return d_data.cols();
}

inline double Data::normalization() const
{
  return d_normalization;
}
