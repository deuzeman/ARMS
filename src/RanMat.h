#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

#include <Data.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Point.h>
#include <Random/stocc.h>

class RanMat
{
  public:
    typedef Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > MCD;
    typedef Eigen::Array< double, Eigen::Dynamic, Eigen::Dynamic> ADD;
    typedef Eigen::Array< unsigned, Eigen::Dynamic, Eigen::Dynamic> UDD;
  
  private:
    size_t const d_N;
    size_t const d_nu;
    double const d_scale;
    
    size_t const d_eigMin;
    size_t const d_numEigs;

    int d_rank;
    int d_nodes;
    size_t d_samples;
    StochasticLib1 d_rstream;
    
    MCD d_Z;
    MCD d_M;

    MCD d_A;
    MCD d_B;
    MCD d_W;
    
    MCD d_gamma_5;

    Eigen::SelfAdjointEigenSolver< MCD > d_slv;
    
  public:
    RanMat(size_t const N, size_t const nu, int const eigMin, int const eigMax);
    ~RanMat();

    void calculate(Point const &params, size_t iter);
    size_t *discretize(double const *breaks, int eigMin, size_t const levels, size_t const eigs) const;

    size_t eigToIndex(int eig) const;
    
    size_t const &numEigs() const;
    size_t const &eigMin() const;
    size_t const &numSamples() const;
    int const &nodes() const;
    int const &rank() const;
    size_t const &N() const;
    size_t const &nu() const;
    
    double const *result() const;
    double *d_result;
};

inline int const &RanMat::nodes() const
{
  return d_nodes;
}

inline int const &RanMat::rank() const
{
  return d_rank;
}

inline double const *RanMat::result() const
{
  return d_result;
}

inline size_t RanMat::eigToIndex(int eig) const
{
  return (eig > 0) ? (d_N + eig - 1) : (d_N + eig);
}

inline size_t const &RanMat::numEigs() const
{
  return d_numEigs;
}

inline size_t const &RanMat::numSamples() const
{
  return d_samples;
}

inline size_t const &RanMat::eigMin() const
{
  return d_eigMin;
}

inline size_t const &RanMat::N() const
{
  return d_N;
}

inline size_t const &RanMat::nu() const
{
  return d_nu;
}
