#pragma once

#include <mpi.h>
#include <Simplex.h>

class Minim
{
  double d_alpha = 1.0;
  double d_gamma = 2.0;
  double d_rho   = -0.5;
  double d_sigma = 0.5;
  
  Simplex      &d_simplex;
  
  public:
    Minim(Simplex &simplex);
    
    Simplex const &reduce();
    double evalPoint(size_t idx);
};

inline Minim::Minim(Simplex &simplex)
  : d_alpha(1.0), d_gamma(2.0), d_rho(-0.5), d_sigma(0.5), d_simplex(simplex)
{}
