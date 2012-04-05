#pragma once

#include <Eigen/Eigen.h>
#include <RanMat/RanMat.h>

class Minim
{
  Data const   *d_data;
  size_t const  d_N;
  int const     d_nu;
  bool const    d_kol;
  
  RanMat        d_engine;
  
  ArrayXXd d_simplex;
  ArrayXd  d_values;
   
  public:
    Minim(Data const *data, size_t const N, int const nu, bool const kol);
    
    void reduce(d_simplex);
};



inline Minim::Minim(Data const *data, size_t const N, int const nu, bool const kol)
: d_data(data), d_N(N), d_nu(nu), d_kol(kol), d_calc()
{
  
}
