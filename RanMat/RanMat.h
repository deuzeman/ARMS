#pragma once

#include <complex>
#include <cstdlib>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

typedef Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > MCD;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ADD;

struct params
{
  size_t N;
  size_t nu;

  double m;
  double a6;
  double a7;
  double a8;

  size_t iter;
  std::string output;

  params();
  params(double const m, double const a6, double const a7, double const a8);
};

inline params::params()
    : N(20),
    nu(0),
    m(0.05),
    a6(0.0),
    a7(0.0),
    a8(0.0),
    iter(10000),
    output("arms.out")
{}

inline params::params(double const m, double const a6, double const a7, double const a8)
    : N(20),
    nu(0),
    m(m),
    a6(a6),
    a7(a7),
    a8(a8),
    iter(10000),
    output("arms.out")
{}

params *parseInput(char const *filename);

class RanMat
{
    size_t const d_N;
    size_t const d_nu;
    size_t const d_nEigs;

    double d_m;
    double d_a6;
    double d_a7;
    double d_a8;

    double const d_scale;

    MCD d_Z;
    MCD d_M;

    MCD d_A;
    MCD d_B;
    MCD d_W;

    ADD d_result;

    Eigen::SelfAdjointEigenSolver< MCD > d_slv;
    StochasticLib1 d_rstream;

  public:
    RanMat(params const *par, int const seed);
    void changeParams(Eigen::VectorXd const &par);

    void calculate(size_t const iter, bool const extend = false);
    double const &result(size_t const nSam, size_t const nEig) const;

    double avResult(size_t const nEig) const;
    double avRatio(size_t const num, size_t const den) const;
};

inline double const &RanMat::result(size_t const nSam, size_t const nEig) const
{
  return d_result.coeffRef(nSam, nEig);
}

inline double RanMat::avResult(size_t const nEig) const
{
  return d_result.col(nEig).mean();
}

inline double RanMat::avRatio(size_t const num, size_t const den) const
{
  return d_result.col(num).cwiseQuotient(d_result.col(den)).mean();
}