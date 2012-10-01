#pragma once

#include <mpi.h>
#include <Simplex.h>

class Minim
{
  double d_alpha;
  double d_gamma;
  double d_rho;
  double d_sigma;
  
  Simplex d_simplex;
  
  public:
    Minim(Data &data, Params &params);
    
    Simplex const &reduce();
    double evalPoint(size_t idx);
};

// NOTE Test: start out using the average, rather than the KS D value.

inline Minim::Minim(Data &data, Params &params)
  : d_alpha(1.0), d_gamma(2.0), d_rho(-0.5), d_sigma(0.5), d_simplex(data, params, Comparator::AVE)
{}
