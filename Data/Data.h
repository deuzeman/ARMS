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

  Eigen::ArrayXXd d_data;

  current mutable d_cur;

  void bootstrap(size_t const nBoot, Eigen::ArrayXXd const &data) const;

  Eigen::ArrayXXd mutable d_summary;

  public:
    Data(char const *filename);

    Eigen::ArrayXXd const &average(size_t const nBoot) const;
    Eigen::ArrayXXd const &ratios(size_t const nBoot) const;
    Eigen::ArrayXXd const &cumulant() const;

    int minEv() const;
    int maxEv() const;
};

inline int Data::minEv() const
{
  return d_minEv;
}

inline int Data::maxEv() const
{
  return d_maxEv;
}
