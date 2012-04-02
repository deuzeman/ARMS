#include <RanMat.h>

#include <algorithm>
#include <iostream>

RanMat::RanMat(size_t const N, size_t const nu, int const nEig_min, int const nEig_max, size_t const nDet, size_t const seed, bool const bootstrap)
  : d_N(N),
    d_nu(nu),
    d_nEig_min(nEig_min),
    d_nEig_max(nEig_max),
    d_nDet(nDet),
    d_seed(seed),
    d_scale(1.0 / std::sqrt(static_cast< double >(2 * d_N + d_nu))),
    d_Z(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_gamma_5(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_A(d_N + d_nu, d_N + d_nu),
    d_B(d_N, d_N),
    d_W(d_N + d_nu, d_N),
    d_slv(2 * d_N + d_nu),
    d_slv_cmplx(2 * d_N + d_nu),
    d_result(2 * d_N + d_nu, d_nEig_max - d_nEig_min),
    d_rstream(d_seed),
    d_bootstrap(bootstrap)
{
  for (int x = 0; x < d_N + d_nu; ++x)
    d_gamma_5(x, x) = 1.0;
  for (int x = d_N + d_nu; x < (2 * d_N + d_nu); ++x)
    d_gamma_5(x, x) = -1.0;
}

void RanMat::calculate(Eigen::ArrayXd const &params, size_t const iter, bool const extend)
{
  size_t offset = 0;
  double const static sqrt1_2 = std::sqrt(2) / 2.0;
  double const static sqrt8   = std::sqrt(8);
  double const scale2  = d_scale * d_scale;
  double const m  = params.coeffRef(0) * scale2;
  double const a6 = params.coeffRef(1) * scale2 * sqrt8;
  double const a7 = params.coeffRef(2) * scale2 * sqrt8;
  double const a8 = params.coeffRef(3) * scale2;
  double const sigma = params.coeffRef(4);

  if (extend)
  {
    offset = d_result.rows();
    d_result.conservativeResize(offset + iter, d_nEig_max - d_nEig_min);
    d_det.conservativeResize(offset + iter);
  }
  else
  {
    d_result.resize(iter, d_nEig_max - d_nEig_min);
    d_det.resize(iter);
  }

  for (size_t ctr = offset; ctr < offset + iter; ++ctr)
  {
    // Initialize block A
    for (int x = 0; x < d_N + d_nu; ++x)
      for (int y = 0; y < d_N + d_nu; ++y)
        d_A(x, y) = std::complex< double >(d_rstream.Normal(0.0, a8), d_rstream.Normal(0.0, a8));
    d_A = sqrt1_2 * (d_A + d_A.adjoint()).eval();

    // Initialize block B
    for (int x = 0; x < d_N; ++x)
      for (int y = 0; y < d_N; ++y)
        d_B(x, y) = std::complex< double >(d_rstream.Normal(0.0, a8), d_rstream.Normal(0.0, a8));
    d_B = sqrt1_2 * (d_B + d_B.adjoint()).eval();

    for (int x = 0; x < d_N + d_nu; ++x)
      for (int y = 0; y < d_N; ++y)
        d_W(x, y) = std::complex< double >(d_rstream.Normal(0.0, d_scale), d_rstream.Normal(0.0, d_scale));

    d_Z << d_A, d_W, d_W.adjoint(), d_B;

    // The following part now adds the d_W6 and d_W7 terms, as well as the mass term
    double mfac = d_rstream.Normal(m, a6);
    double ifac = d_rstream.Normal(0.0, a7);
    d_Z += mfac * d_gamma_5 + ifac * MCD::Identity(2 * d_N + d_nu, 2 * d_N + d_nu);
    d_slv.compute(d_Z, Eigen::EigenvaluesOnly);
    int shift = 0;
    while (d_slv.eigenvalues()[d_N + shift - 1] >= 0)
      --shift;
    while (d_slv.eigenvalues()[d_N + shift] < 0)
      ++shift;
    d_result.row(ctr) = ((2 * d_N + d_nu) / sigma) * d_slv.eigenvalues().segment(static_cast< int >(d_N + shift) + d_nEig_min, d_nEig_max - d_nEig_min);
    d_det[ctr] = (d_nDet > 0) ? (((2 * d_N + d_nu) / sigma) * d_slv.eigenvalues().segment(d_N -(d_nDet / 2), d_nDet)).prod() : 0;
  }
}

void RanMat::calcreal(Eigen::ArrayXd const &params)
{
  double const static sqrt1_2 = std::sqrt(2) / 2.0;
  double const static sqrt8   = std::sqrt(8);
  double const scale2  = d_scale * d_scale;
  double const m  = params.coeffRef(0) * scale2;
  double const a6 = params.coeffRef(1) * scale2 * sqrt8;
  double const a7 = params.coeffRef(2) * scale2 * sqrt8;
  double const a8 = params.coeffRef(3) * scale2;

  // Initialize block A
  for (int x = 0; x < d_N + d_nu; ++x)
    for (int y = 0; y < d_N + d_nu; ++y)
      d_A(x, y) = std::complex< double >(d_rstream.Normal(0.0, a8), d_rstream.Normal(0.0, a8));
  d_A = sqrt1_2 * (d_A + d_A.adjoint()).eval();

  // Initialize block B
  for (int x = 0; x < d_N; ++x)
    for (int y = 0; y < d_N; ++y)
      d_B(x, y) = std::complex< double >(d_rstream.Normal(0.0, a8), d_rstream.Normal(0.0, a8));
  d_B = sqrt1_2 * (d_B + d_B.adjoint()).eval();

  for (int x = 0; x < d_N + d_nu; ++x)
    for (int y = 0; y < d_N; ++y)
      d_W(x, y) = std::complex< double >(d_rstream.Normal(0.0, d_scale), d_rstream.Normal(0.0, d_scale));
  d_Z << d_A, d_W, -1 * d_W.adjoint(), d_B;

  // The following part now adds the d_W6 and d_W7 terms, as well as the mass term
  double mfac = d_rstream.Normal(m, a6);
  double ifac = d_rstream.Normal(0.0, a7);
  d_Z += mfac * MCD::Identity(2 * d_N + d_nu, 2 * d_N + d_nu) + ifac * d_gamma_5;
  d_slv_cmplx.compute(d_Z, true);
}

Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::bootAver(size_t const nBoot) const
{
  d_average.resize(2, d_result.cols());
  d_average.row(0) = d_result.colwise().mean();

  if (d_bootstrap)
  {
    Eigen::ArrayXXd bootSamp(d_result.rows(), d_result.cols());
    Eigen::ArrayXXd bootHist(nBoot, d_result.cols());
    CRandomSFMT sampler(time(0));

    for (size_t bootCtr = 0; bootCtr < nBoot; ++bootCtr)
    {
      for (size_t sampCtr = 0; sampCtr < d_result.rows(); ++sampCtr)
        bootSamp.row(sampCtr) = d_result.row(sampler.IRandomX(0, d_result.rows() - 1));
      bootHist.row(bootCtr) = (bootSamp.colwise().mean());
    }

    // Calculate the standard deviation from this.
    d_average.row(1) = bootHist.square().colwise().mean();
    d_average.row(1) -= bootHist.colwise().mean().square();
    d_average.row(1) = std::sqrt(d_average.row(1));
  }

  return d_average;
}

Eigen::Array< double, 2, Eigen::Dynamic > const &RanMat::bootRat(size_t const nBoot) const
{
  d_ratios.resize(2, d_result.cols() * (d_result.cols() - 1) / 2);
  Eigen::ArrayXXd ratFull(d_result.rows(), d_result.cols() * (d_result.cols() - 1) / 2);

  size_t ctr = 0;
  for (size_t num = 0; num < d_result.cols() - 1; ++num)
    for (size_t den = num + 1; den < d_result.cols(); ++den)
      ratFull.col(ctr++) = (d_result.col(num) / d_result.col(den));

  d_ratios.row(0) = ratFull.colwise().mean();

  if (d_bootstrap)
  {
    Eigen::ArrayXXd bootSamp(ratFull.rows(), ratFull.cols());
    Eigen::ArrayXXd bootHist(nBoot, ratFull.cols());
    CRandomSFMT sampler(time(0));

    for (size_t bootCtr = 0; bootCtr < nBoot; ++bootCtr)
    {
      for (size_t sampCtr = 0; sampCtr < ratFull.rows(); ++sampCtr)
	bootSamp.row(sampCtr) = ratFull.row(sampler.IRandomX(0, ratFull.rows() - 1));
      bootHist.row(bootCtr) = (bootSamp.colwise().mean());
    }

    // Calculate the standard deviation from this.
    d_ratios.row(1) = bootHist.square().colwise().mean();
    d_ratios.row(1) -= bootHist.colwise().mean().square();
    d_ratios.row(1) = std::sqrt(d_ratios.row(1));
  }

  return d_ratios;
}

double *RanMat::flat() const
{
  double *pdata = new double[d_result.rows() *  d_result.cols()];
  for (size_t row = 0; row <  d_result.rows(); ++row)
    for (size_t col = 0; col <  d_result.cols(); ++col)
      pdata[row *  d_result.cols() + col] =  d_result.coeff(row, col);
  std::sort(pdata, pdata +  d_result.rows() *  d_result.cols());
  return pdata;
}

double *RanMat::flat(size_t const &col) const
{
  double *pdata = new double[d_result.rows()];
  for (size_t row = 0; row <  d_result.rows(); ++row)
    pdata[row] = d_result.coeff(row, col);
  std::sort(pdata, pdata +  d_result.rows());
  return pdata;
}
