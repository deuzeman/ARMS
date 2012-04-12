#pragma once

#include <RanMat.h>

class Discretizer
{
  size_t d_numEigs;
  size_t d_levels;
  size_t d_numBlocks;
  double const *d_breaks;
  
  size_t *d_data;
  double *d_hist;
  double *d_acc;
  
  double *d_cum;
  double *d_block;
  
  size_t d_sampTotal;
  
  public:
    Discretizer(double const *breaks, size_t numBreaks, size_t numEigs, size_t blocks);
    ~Discretizer();
    
    void clear();
    void calculate(RanMat const &ranmat, double const *breaks, size_t numBreaks, size_t blocks);
    double const &operator()(size_t const &eig, size_t const &level) const;
    inline double const &operator()(size_t const &eig, size_t const &level, size_t const &block) const;
};

inline double const &Discretizer::operator()(size_t const &eig, size_t const &level) const
{
  return d_cum[eig * d_levels + level];
}

inline double const &Discretizer::operator()(size_t const &eig, size_t const &level, size_t const &block) const
{
  return d_block[(block * d_numEigs + eig) * d_levels + level];
}
