#pragma once

#include <numeric>
#include <RanMat.h>

class Discretizer
{
  int    d_nodes;
  
  size_t d_numEigs;
  size_t d_levels;
  size_t d_numBlocks;
  size_t d_numBreaks;
  double const *d_breaks;
  
  size_t *d_data;
  double *d_hist;
  double *d_acc;
  
  double *d_cum;
  double *d_block;
  double d_blockFac;
  
  size_t d_sampTotal;
  
  double *d_avBlocks;
  
  public:
    Discretizer(double const *breaks, size_t numBreaks, size_t numEigs, size_t blocks);
    ~Discretizer();
    
    void clear();
    void calculate(RanMat const &ranmat);
    double const &operator()(size_t const &eig, size_t const &level) const;
    double operator()(size_t const &eig, size_t const &level, size_t const &block) const;
    
    double average(size_t const &eig) const;
    double average(size_t const &eig, size_t const &block) const;
};

inline double const &Discretizer::operator()(size_t const &eig, size_t const &level) const
{
  return d_cum[eig * d_levels + level];
}

inline double Discretizer::operator()(size_t const &eig, size_t const &level, size_t const &block) const
{
  return d_blockFac * (d_cum[eig * d_levels + level] - d_block[(block * d_numEigs + eig) * d_levels + level]);
}

inline double Discretizer::average(size_t const &eig) const
{
  double res = 0.0;
  for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
    res += d_avBlocks[eig * d_numBlocks + bIdx];
  return res /= d_numBlocks;
}

inline double Discretizer::average(size_t const &eig, size_t const &block) const
{
  double res = 0.0;
  for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
  {
    if (bIdx == block)
      continue;
    res += d_avBlocks[eig * d_numBlocks + bIdx];
  }
  return res /= (d_numBlocks - 1);
}