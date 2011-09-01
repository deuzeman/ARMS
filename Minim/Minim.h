#pragma once

#include <cmath>
#include <iostream>

#include <Data.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <RanMat.h>

// Uses a combination of Powell's method with Brent to perform derivative-free minimization

class Minim
{
  int d_nEig_min;
  int d_nEig_max;
  RanMat *d_generator;
  RanMat *d_secondary;
  Data const &d_data;

  Eigen::ArrayXd d_res;

  void switchGen();
  double round(double const value, double const tol);
  
  template< typename MT >
  MT &round(MT &vec, double const tol);
  
  double phi(double const start, double const edge, size_t const iter, Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir,
             double const value, double const error, double const tol, double * const bval, double * const berr);

  public:
    Minim(Data const &data, size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const seed);
    ~Minim();

    double chiSq(Eigen::ArrayXd const &pars, size_t const iter, double * const error = 0, size_t const nBoot = 1000, bool const extend = false);
    double brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, size_t const rmIters = 50000, double const tol = 1e-8);
    void powell(Eigen::ArrayXd &start, Eigen::ArrayXXd &bounds, size_t rmIters = 50000, size_t const powIters = 100, double const tol = 1e-8);
};

inline Minim::Minim(Data const &data, size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const seed)
  : d_nEig_min(nEig_min), d_nEig_max(nEig_max), d_generator(new RanMat(N, nu, nEig_min, nEig_max, 0, seed)), d_secondary(new RanMat(N, nu, nEig_min, nEig_max, 0, seed + 235)), d_data(data)
{}

inline Minim::~Minim()
{
  delete d_generator;
  delete d_secondary;
}

inline void Minim::switchGen()
{
  RanMat *tmp = d_generator;
  d_generator = d_secondary;
  d_secondary = tmp;
}

inline double Minim::round(double const value, double const tol)
{
  return tol * std::floor(0.5 + value / tol);
}

template< typename MT >
MT &Minim::round(MT &vec, double const tol)
{
  vec.array() /= tol;
  vec.array() += 0.5;
  for (size_t row = 0; row < vec.rows(); ++row)
    for (size_t col = 0; col < vec.cols(); ++col)
      vec(row, col) = std::floor(vec(row, col));
  vec.array() *= tol;
  return vec;
}
