#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>

// Uses a combination of Powell's method with Brent to perform derivative-free minimization

class Minim
{
  public:
    Eigen::MatrixXd brackets;

    Eigen::ArrayXd  ratData;
    Eigen::ArrayXd  sdRatData;
    
    size_t iter;
    size_t maxBrIter;
    
    double tol;
    
  private:
    Eigen::MatrixXd d_dirs;
    Eigen::MatrixXd d_pars;
    Eigen::ArrayXd  d_vals;

    RanMat d_generator;

  public:
    Minim(size_t const dim, size_t const nEigs = 12);

  private:
    void brent(size_t const idx);
    void predict(Eigen::VectorXd const &pars, size_t const iter, Eigen::ArrayXd &eval);
};

inline Minim::Minim(size_t const dim, size_t const nEigs = 12)
    : brackets(Eigen::MatrixXd::Zeros(dim, 2))
    d_dirs(Eigen::MatrixXd::Identity(dim)),
    d_pars(Eigen::MatrixXd::Zeros(dim, dim + 1)),
    d_eval(Eigen::MatrixXd::Zeros(nEigs * (nEigs - 1) / 2))
{}


