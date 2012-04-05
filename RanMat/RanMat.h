#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
    StochasticLib1 d_rstream;
    
    MCD d_Z;
    MCD d_M;

    MCD d_A;
    MCD d_B;
    MCD d_W;
    
    MCD d_gamma5;

    Eigen::SelfAdjointEigenSolver< MCD > d_slv;
    
    MCD d_result;
    UDD d_result_discrete;

  public:
    RanMat(int const eigMin, int const eigMax);

    void calculate(Eigen::ArrayXd const &params, size_t iter, bool const extend = false);
    void discretize(UDD const &breaks);
    
    MCD const &result() const;
    UDD const &result_discrete() const;
    size_t const &numEigs() const;
    size_t const nodes() const;
};

void histogram(RanMat::ADD *result, RanMat const &sim, double min, double max, double bins);

inline RanMat::MCD const &RanMat::result() const
{
  return d_result;
}

inline size_t const &RanMat::numEigs() const
{
  return d_numEigs;
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
