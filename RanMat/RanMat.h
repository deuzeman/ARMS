#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
    
    MCD d_gamma5;

    Eigen::SelfAdjointEigenSolver< MCD > d_slv;
    
    double *d_result;
    double *d_resultDiscrete;
    bool d_isDiscretized;

  public:
    RanMat(size_t const N, size_t const nu, int const eigMin, int const eigMax);
    ~RanMat();

    void calculate(Point const &params, size_t iter, bool const extend = false);
    void discretizedouble const *breaks, int eigMin, size_t const levels, size_t const eigs);

    size_t eigToIndex(int eig) const;
    
    double const *result() const;
    size_t const *result_discrete() const;
    size_t const &numEigs() const;
    size_t const &numSamples() const;
    size_t const &nodes() const;
    size_t const &N() const;
    int const &nu() const;
    bool const &isDiscretized() const;
};

void kolmogorov(StatVal *kol, RanMat const &sim, Data const &data);

inline RanMat::MCD const &RanMat::result() const
{
  return d_result;
}

inline double const *RanMat::resultDiscrete() const
{
  return d_resultDiscrete;
}

inline size_t RanMat::eigToIndex(int eig) const
{
  return (eig > 0) ? (d_N + eig - 1) : (d_N + eig);
}

inline bool const &isDiscretized() const
{
  return d_is_discretized;
}

inline size_t const &RanMat::numEigs() const
{
  return d_numEigs;
}

inline size_t const &RanMat::numSamples() const
{
  return d_result.rows();
}

inline size_t const &RanMat::eigMin() const
{
  return d_eigMin;
}

inline size_t const &RanMat::nodes() const
{
  return d_nodes;
}

inline size_t const &RanMat::N() const
{
  return d_N;
}

inline int const &RanMat::nu() const
{
  return d_nu;
}
