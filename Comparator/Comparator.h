#pragma once

#include <mpi.h>

#include <Data.h>
#include <Point.h>
#include <RanMat.h>


class Comparator
{
  RanMat &d_ranmat;
  
  double *d_breaks;
  size_t  d_levels;
  double  d_inc;
  size_t  d_eigs;
  int     d_minEv;
 
  double  d_relprec
  size_t  d_blocks;
  
  int     d_rank;
  int     d_nodes;
  
  public:
    Comparator(Data &data, RanMat &ranmat);
    void setPrecision(double relprec);
    void setBlocks(size_t blocks);
    
    double kolmogorov(Point const &point);
    
    size_t roundToBlocks(size_t in) const;
};

inline Comparator::Comparator(Data &data, RanMat &ranmat)
: d_ranmat(ranmat), 
  d_breaks(data.flatPerColumn()), 
  d_levels(data.numSamples()), 
  d_inc(1.0 / d_levels), 
  d_eigs(data.numCols()), 
  d_minEv(data.minEv()),
  d_relprec(1.0e-4), 
  d_blocks(50)
{
  MPI_Get_rank(&d_rank, MPI_COMM_WORLD);
  MPI_Comm_size(&d_nodes, MPI_COMM_WORLD);
}


inline void Comparator::setPrecision(double relprec)
{
  d_relprec = relprec;
}

inline void Comparator::setBlocks(size_t blocks)
{
  d_blocks = blocks;
}

inline size_t Comparator::roundToBlocks(size_t in) const
{
  // Give an alignment such that we always have a multiple of blocks samples per node.
  return ((((in - 1) / (d_blocks * d_nodes)) + 1) * (d_blocks * d_nodes));
}