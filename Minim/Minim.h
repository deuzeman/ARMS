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
  RanMat d_generator;
  Data const &d_data;

  Eigen::ArrayXd d_res;

  public:
    Minim(Data const &data, size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const seed);

    double chiSq(Eigen::ArrayXd const &pars, size_t const iter);
    double brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, size_t const rmIters = 50000, double const tol = 1e-8);
    void powell(Eigen::VectorXd const &start, Eigen::ArrayXXd const &bounds, size_t const rmIters = 50000, size_t const powIters = 100, double const tol = 1e-8);
};

inline Minim::Minim(Data const &data, size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const seed)
  : d_nEig_min(nEig_min), d_nEig_max(nEig_max), d_generator(N, nu, nEig_min, nEig_max, 0, seed), d_data(data)
{}
