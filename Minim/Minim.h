#pragma once

#include <mpi.h>
#include <Eigen/Eigen.h>
#include <Params/Params.h>
#include <Point/Point.h>
#include <RanMat/RanMat.h>

struct Simplex
{
  Point   points[6];
  StatVal values[6];
  
  Simplex(Point const &p);
};

inline Simplex::Simplex(Point const &p)
{
  int rank;
  if (rank == 0)
  {
    seed(time());
    points[0] = p;
    // Plug in random shifts for the other values
    for (size_t idx = 1; idx < 6; ++idx)
    {
      points[idx].sigma = p.sigma * (1 + (0.5 * ((rand() / RANDMAX) - 1)));
      points[idx].m     = p.m     * (1 + (0.5 * ((rand() / RANDMAX) - 1)));
      points[idx].a6    = (std::abs(a6) > 1e-8) ? p.a6 * (1 + (0.5 * ((rand() / RANDMAX) - 1))) : 0.0;
      points[idx].a7    = (std::abs(a7) > 1e-8) ? p.a7 * (1 + (0.5 * ((rand() / RANDMAX) - 1))) : 0.0;
      points[idx].a8    = (std::abs(a8) > 1e-8) ? p.a8 * (1 + (0.5 * ((rand() / RANDMAX) - 1))) : 0.0;
    }
  }
  MPI_Bcast(points, 30, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

class Minim
{
  Data const   *d_data;
  size_t const  d_N;
  int const     d_nu;
  
  Simplex       d_simplex;
  RanMat        d_engine;
  
  public:
    Minim(Data const *data, FitParams const *params);
    
    Point *reduce(Point *initial, double prec);
    double evalPoint(size_t idx) const;
    void   sortPoints(size_t *ranking, Simplex *simplex);
    void   sortLocation(size_t idx, Simplex *simplex);
};

inline Minim::Minim(Data const *data, FitParams const *params)
: d_data(data), d_N(N), d_nu(nu), d_simplex(params->p), d_engine(params->N, params->nu, data->eigMin, data->eigMax)
{}
