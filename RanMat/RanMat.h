#pragma once

#include <complex>
#include <cstdlib>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

typedef Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > MCD;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ADD;

class RanMat
{
  size_t const d_N;
  size_t const d_nu;
  size_t const d_nEig_min;
  size_t const d_nEig_max;
  size_t const d_nDet;
  size_t const d_seed;
  bool const d_bootstrap;

  double const d_scale;

  MCD d_Z;
  MCD d_M;

  MCD d_A;
  MCD d_B;
  MCD d_W;

  ADD d_result;

  Eigen::ArrayXd d_det;
  Eigen::Array< double, 2, Eigen::Dynamic > mutable d_ratios;
  Eigen::Array< double, 2, Eigen::Dynamic > mutable d_average;

  Eigen::SelfAdjointEigenSolver< MCD > d_slv;
  StochasticLib1 d_rstream;

  void RanMat::bootAver(size_t const nBoot) const;
  void RanMat::bootRat(size_t const nBoot) const;

  public:
    RanMat(size_t const N, size_t const nu, size_t const nEig_min, size_t const nEig_max, size_t const nDet, size_t const seed, bool const bootstrap = false);

    void calculate(Eigen::ArrayXd const &params, size_t const iter, bool const extend = false);
    double const &result(size_t const nSam, size_t const nEig) const;
    double const &det(size_t const nSam) const;

    Eigen::Array< double, 2, Eigen::Dynamic > const &average() const;
    Eigen::Array< double, 2, Eigen::Dynamic > const &ratios() const;

    Eigen::Array< double, Eigen::Dynamic, Eigen::Dynamic > const &samples() const;
};

inline double const &RanMat::result(size_t const nSam, size_t const nEig) const
{
  return d_result.coeffRef(nSam, nEig);
}

inline double const &RanMat::det(size_t const nSam) const
{
  return d_det.coeffRef(nSam);
}

inline Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::average() const
{
  return d_average;
}

inline Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::ratios() const
{
  return d_ratios;
}

inline Eigen::Array< double, Eigen::Dynamic, Eigen::Dynamic > const &RanMat::samples() const
{
  return d_result;
}
