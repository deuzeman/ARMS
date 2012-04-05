#pragma once

#include <cstdlib>
#include <fstream>
#include <sstream>

#include <Eigen/Dense>
#include <Random/sfmt.h>


class Data
{
  enum current
  {
    NONE,
    AVE,
    RAT,
    CUM
  };

  int d_minEv;
  int d_maxEv;
  size_t d_nSamp;
  int d_bootSeed;

  Eigen::ArrayXXd d_data;
  double d_normalization;

  current mutable d_cur;

  void bootstrap(size_t const nBoot, Eigen::ArrayXXd const &data) const;

  Eigen::ArrayXXd mutable d_summary;

  public:
    Data(char const *filename, int bootSeed = 0, bool const parallel = false);

    Eigen::ArrayXXd const &average(size_t const nBoot) const;
    Eigen::ArrayXXd const &cumulant() const;

    int minEv() const;
    int maxEv() const;
    size_t numSamples() const;
    size_t numCols() const;
    int bootSeed() const;

    double normalization() const;

    double *flat() const;
    double *flat(size_t const &col) const;
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

inline int Data::bootSeed() const
{
  return d_bootSeed;
}

inline double Data::normalization() const
{
  return d_normalization;
}
