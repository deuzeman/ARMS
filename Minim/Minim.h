#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <mpi.h>

#include <Data.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <RanMat.h>

// Uses a combination of Powell's method with Brent to perform derivative-free minimization

struct Control
{
  int order[2];
  double params[5];
};

class Minim
{
  // Code for the command and control structure
  MPI_Datatype MPI_CONTROL;
  int d_mpirank;
  int d_numnodes;

  size_t d_kol;

  Control d_control;
  unsigned long int *d_countBuffer;
  unsigned long int *d_countBackup;

  RanMat *d_generator;
  Eigen::ArrayXXi *d_hash_gen;

  RanMat *d_secondary;
  Eigen::ArrayXXi *d_hash_sec;

  CRandomSFMT *d_sampler;
  Data const &d_data;

  Eigen::ArrayXd d_res;

  std::ofstream d_log;

  void switchGen();
  void condense();
  void buildHash(bool const extend);
  void clearCountBuffer();
  void addSampleToCount(size_t const idx);
  double round(double const value, double const tol);

  template< typename MT >
  MT &round(MT &vec, double const tol);

  double brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, 
               size_t const rmIters = 1000, size_t const maxIters = 50000, double const tol = 1e-8);
  double phi(double const start, double const edge, size_t const iter, Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir,
             double const value, double const error, double const tol, double * const bval, double * const berr);

  public:
    Minim(Data const &data, size_t const N, size_t const nu, size_t const kol);
    ~Minim();

    void listen();

    double chiSq(Eigen::ArrayXd const &pars, size_t const iter, double * const error = 0, size_t const nBoot = 1000, bool const extend = false);
    void powell(Eigen::ArrayXd &start, Eigen::ArrayXXd &bounds, size_t rmIters = 1000, size_t const maxIters = 50000, size_t const powIters = 100, double const tol = 1e-8);

    std::ofstream &log();
};

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

inline std::ofstream &Minim::log()
{
  return d_log;
}

inline void Minim::clearCountBuffer()
{
  std::fill_n(d_countBuffer, d_data.numCols() * d_data.numSamples(), 0);
}

inline void Minim::addSampleToCount(size_t const idx)
{
  for (size_t col = 0; col < d_data.numCols(); ++col)
    ++d_countBuffer[(*d_hash_gen)(idx, col) + col * d_data.numSamples()];
}