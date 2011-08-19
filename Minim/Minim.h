#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>

// Uses a combination of Powell's method with Brent to perform derivative-free minimization

class Minim
{
  public:
    Eigen::MatrixXd brackets;

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

inline Minim::Minim(size_t const dim)
  : brackets(Eigen::MatrixXd::Zeros(dim, 2))
    d_dirs(Eigen::MatrixXd::Identity(dim)),
    d_pars(Eigen::MatrixXd::Zeros(dim, dim + 1)),
    d_eval(Eigen::MatrixXd::Zeros(nEigs * (nEigs - 1) / 2);
{}

void Minim::brent(size_t const idx)
{
  // Performs Brent's line search in the direction d_dirs(idx) with start d_pars(idx)

  // Determine brackets on the scaling factor...
  Eigen::MatrixXd limits(dim, 2);
  limits.row(0) = (brackets.row(0) - d_pars.col(idx)) / d_dirs.col(idx);
  limits.row(1) = (brackets.row(1) - d_pars.col(idx)) / d_dirs.col(idx); // Can be better, perhaps.

  double min = limits.coeff(0, limits.row(0).minCoeff());
  double max = limits.coeff(1, limits.row(1).maxCoeff());
}

void Minim::chiSq(Eigen::VectorXd const &pars, size_t const iter, Eigen::ArrayXd const &ratData, Eigen::ArrayXd const &sdRatData)
{
  static Eigen::ArrayXd eval;

  d_generator.changeParams(params(pars));
  d_generator.calculate(iter);

  size_t resCtr = 0;
  for (size_t num = 0; num < nEigs - 1; ++num)
    for (size_t den = num; den < nEigs; ++den)
      eval.coeff(resCtr++) = d_generator.avRatio(num, den);

  (eval - ratData)
}
