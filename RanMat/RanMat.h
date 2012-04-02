#pragma once

#include <complex>
#include <cstdlib>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Random/stocc.h>

typedef Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > MCD;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ADD;
typedef Eigen::Array< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic> ACD;

class RanMat
{
  size_t const d_N;
  size_t const d_nu;
  int const d_nEig_min;
  int const d_nEig_max;
  size_t const d_nDet;
  size_t const d_seed;
  bool const d_bootstrap;

  double const d_scale;

  MCD d_Z;
  MCD d_M;
  MCD d_gamma_5;

  MCD d_A;
  MCD d_B;
  MCD d_W;

  ADD d_result;
  ACD d_result_cmplx;

  Eigen::ArrayXd d_det;
  Eigen::Array< double, 2, Eigen::Dynamic > mutable d_ratios;
  Eigen::Array< double, 2, Eigen::Dynamic > mutable d_average;

  Eigen::SelfAdjointEigenSolver< MCD > d_slv;
  Eigen::ComplexEigenSolver< MCD > d_slv_cmplx;
  StochasticLib1 d_rstream;

  Eigen::Array< double, 2, Eigen::Dynamic > const &bootAver(size_t const nBoot) const;
  Eigen::Array< double, 2, Eigen::Dynamic > const &bootRat(size_t const nBoot) const;

  public:
    RanMat(size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const nDet, size_t const seed, bool const bootstrap = false);

    void calculate(Eigen::ArrayXd const &params, size_t const iter, bool const extend = false);
    void calcreal(Eigen::ArrayXd const &params);
    double const &result(size_t const nSam, size_t const nEig) const;
    std::complex< double > const &result_real(size_t const nSam, size_t const nEig) const;
    double const &det(size_t const nSam) const;

    Eigen::Array< double, 2, Eigen::Dynamic > const &average() const;
    Eigen::Array< double, 2, Eigen::Dynamic > const &ratios() const;

    ADD const &samples() const;
    MCD const &gamma_5() const;
    Eigen::Matrix< std::complex< double >, Eigen::Dynamic, 1 > const &eigenvalues_cmplx() const;
    Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > const &eigenvectors_cmplx() const;

    double *flat() const;
    double *flat(size_t const &col) const;
};

inline double const &RanMat::result(size_t const nSam, size_t const nEig) const
{
  return d_result.coeffRef(nSam, nEig);
}

inline std::complex< double > const &RanMat::result_real(size_t const nSam, size_t const nEig) const
{
  return d_result_cmplx.coeffRef(nSam, nEig);
}

inline double const &RanMat::det(size_t const nSam) const
{
  return d_det.coeffRef(nSam);
}

inline Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::average() const
{
  return bootAver(1000);
}

inline Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::ratios() const
{
  return bootRat(1000);
}

inline ADD const &RanMat::samples() const
{
  return d_result;
}

inline MCD const &RanMat::gamma_5() const
{
  return d_gamma_5;
}

inline Eigen::Matrix< std::complex< double >, Eigen::Dynamic, 1 > const &RanMat::eigenvalues_cmplx() const
{
  return d_slv_cmplx.eigenvalues();
}

inline Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > const &RanMat::eigenvectors_cmplx() const
{
  return d_slv_cmplx.eigenvectors();
}
