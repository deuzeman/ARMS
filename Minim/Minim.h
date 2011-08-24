#pragma once

#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <RanMat.h>

// Uses a combination of Powell's method with Brent to perform derivative-free minimization

class Minim
{
  size_t d_nEigs;
  RanMat d_generator;
  Eigen::ArrayXXd const &d_data;

  Eigen::ArrayXd d_res;

  public:
    Minim(Eigen::ArrayXXd const &data, size_t const N, size_t const nu, size_t const nEigs, size_t const seed);

  private:
    double chiSq(Eigen::ArrayXd const &pars, size_t const iter);
    double brent(Eigen::ArrayXd const &center, Eigen::ArrayXd const &dir, Eigen::ArrayXXd const &bounds, size_t const rmIters = 50000, double const tol = 1e-8);
    void powell(Eigen::VectorXd const &start, Eigen::ArrayXXd const &bounds, size_t const rmIters = 50000, size_t const powIters = 100, double const tol = 1e-8);
};

inline Minim::Minim(Eigen::ArrayXXd const &data, size_t const N, size_t const nu, size_t const nEigs, size_t const seed)
  : d_nEigs(nEigs), d_generator(N, nu, nEigs, 0, seed), d_data(data)
{}
