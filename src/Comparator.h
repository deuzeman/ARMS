#pragma once

#include <mpi.h>

#include <Data.h>
#include <Discretizer.h>
#include <Log.h>
#include <Params.h>
#include <Point.h>
#include <RanMat.h>

class Comparator
{
  RanMat  d_ranmat;
  
  double *d_aver;
  size_t  d_levels;
  double  d_inc;
  size_t  d_eigs;
  int      d_minEv;
 
  double  d_prec;
  size_t  d_blocks;
  
  int     d_rank;
  int     d_nodes;

  // The following provides scratch space
  Discretizer d_disc;
  double *d_jack;  
  double *d_cumdist;
  
  public:
    Comparator(Data &data, Params &params);
    ~Comparator();
    
    double kolmogorov(Point const &point);
    double averages(Point const &point);
    
    size_t roundToBlocks(size_t in) const;
};

inline size_t Comparator::roundToBlocks(size_t in) const
{
  // Give an alignment such that we always have a multiple of blocks samples per node.
  return ((((static_cast< int >(in) - 1) / (d_blocks * d_nodes)) + 1) * (d_blocks * d_nodes));
}
