#pragma once

#include <numeric>
#include <RanMat.h>

class Discretizer
{
  int    d_nodes;
  
  unsigned long d_numEigs;
  unsigned long d_numBlocks;
  unsigned long d_numLevels;
  unsigned long d_sampTotal;

  double const *d_breaks;
  
  double *d_sum_blocks;
  unsigned long *d_hist_blocks;
  
  double *d_mean_blocks;
  double *d_cum_blocks;
  
  public:
    Discretizer(double const *breaks, unsigned long numBreaks, unsigned long numEigs, unsigned long blocks);
    ~Discretizer();
    
    void clear();
    void calculate(RanMat const &ranmat);
    double operator()(unsigned long const &eig, unsigned long const &level) const;
    double operator()(unsigned long const &eig, unsigned long const &level, unsigned long const &block) const;
    
    double average(unsigned long const &eig) const;
    double average(unsigned long const &eig, unsigned long const &block) const;
};

inline double Discretizer::operator()(unsigned long const &eig, unsigned long const &level, unsigned long const &block) const
{
  return d_cum_blocks[(block * d_numEigs + eig) * d_numLevels + level];
}

inline double Discretizer::average(unsigned long const &eig, unsigned long const &block) const
{
  return d_mean_blocks[block * d_numEigs + eig];
}
