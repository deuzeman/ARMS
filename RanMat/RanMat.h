#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

#include <Data/Data.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Point/Point.h>
#include <Random/stocc.h>


struct StatVal
{
  double value;
  double error;
};

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
    
    double *d_result;
    mutable size_t *d_resultDiscrete;
    mutable bool d_isDiscretized;

  public:
    RanMat(size_t const N, size_t const nu, int const eigMin, int const eigMax);
    ~RanMat();

    void calculate(Point const &params, size_t iter, bool const extend = false);
    size_t *discretize(double const *breaks, int eigMin, size_t const levels, size_t const eigs) const;

    size_t eigToIndex(int eig) const;
    
    double const *result() const;
    size_t const *resultDiscrete() const;
    size_t const &numEigs() const;
    size_t const &eigMin() const;
    size_t const &numSamples() const;
    int const &nodes() const;
    size_t const &N() const;
    size_t const &nu() const;
    bool const &isDiscretized() const;
};

void kolmogorov(StatVal *kol, RanMat const &sim, Data const &data);

inline int const &RanMat::nodes() const
{
  return d_nodes;
}

inline double const *RanMat::result() const
{
  return d_result;
}

inline size_t const *RanMat::resultDiscrete() const
{
  return d_resultDiscrete;
}

inline size_t RanMat::eigToIndex(int eig) const
{
  return (eig > 0) ? (d_N + eig - 1) : (d_N + eig);
}

inline bool const &RanMat::isDiscretized() const
{
  return d_isDiscretized;
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
