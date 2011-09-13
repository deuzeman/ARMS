#include <RanMat.h>

#include <algorithm>
#include <iostream>

RanMat::RanMat(size_t const N, size_t const nu, size_t const nEig_min, size_t const nEig_max, size_t const nDet, size_t const seed, bool const bootstrap)
  : d_N(N),
    d_nu(nu),
    d_nEig_min(nEig_min),
    d_nEig_max(nEig_max),
    d_nDet(nDet),
    d_seed(seed),
    d_scale(0.5 / std::sqrt(static_cast< double >(2 * d_N + d_nu))),
    d_Z(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_M(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_A(d_N + d_nu, d_N + d_nu),
    d_B(d_N, d_N),
    d_W(d_N + d_nu, d_N),
    d_slv(2 * d_N + d_nu),
    d_result(2 * d_N + d_nu, d_nEig_max - d_nEig_min),
    d_rstream(d_seed),
    d_bootstrap(bootstrap)
{}

void RanMat::calculate(Eigen::ArrayXd const &params, size_t const iter, bool const extend)
{
  size_t offset = 0;
  double const &m  = params.coeffRef(0);
  double const &a6 = params.coeffRef(1);
  double const a6_b_m = a6 / m;
  double const &a7 = params.coeffRef(2);
  double const &a8 = params.coeffRef(3);
  double const &scale = params.coeffRef(4);

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

  for (int x = 0; x < d_N + d_nu; ++x)
    d_M(x, x) = m;
  for (int x = d_N + d_nu; x < (2 * d_N + d_nu); ++x)
    d_M(x, x) = -m;

  for (size_t ctr = offset; ctr < offset + iter; ++ctr)
  {
    // Initialize block A
    for (int x = 0; x < d_N + d_nu; ++x)
    {
      d_A(x, x) = std::complex< double >(d_rstream.Normal(0.0, d_scale), 0.0);
      for (int y = x + 1; y < d_N + d_nu; ++y)
      {
        d_A(x, y) = std::complex< double >(d_rstream.Normal(0.0, d_scale), d_rstream.Normal(0.0, d_scale));
        d_A(y, x) = std::conj(d_A(x, y));
      }
    }

    // Initialize block B
    for (int x = 0; x < d_N; ++x)
    {
      d_B(x, x) = std::complex< double >(d_rstream.Normal(0.0, d_scale), 0.0);
      for (int y = x + 1; y < d_N; ++y)
      {
        d_B(x, y) = std::complex< double >(d_rstream.Normal(0.0, d_scale), d_rstream.Normal(0.0, d_scale));
        d_B(y, x) = std::conj(d_B(x, y));
      }
    }

    for (int x = 0; x < d_N + d_nu; ++x)
      for (int y = 0; y < d_N; ++y)
        d_W(x, y) = std::complex< double >(d_rstream.Normal(0.0, d_scale), d_rstream.Normal(0.0, d_scale));

    d_Z << a8 * d_A, d_W, d_W.adjoint(), a8 * d_B;

    // The following part now adds the d_W6 and d_W7 terms, as well as the mass term
    double mfac = d_rstream.Normal(1.0, d_scale * a6_b_m);
    double ifac = d_rstream.Normal(0.0, d_scale * a7);
/*    for (int x = 0; x < (2 * d_N + d_nu); ++x)
      d_Z(x,x) += (d_rstream.Normal(1.0, d_scale * a6_b_m)) * d_M(x, x) + d_rstream.Normal(0.0, d_scale * a7);*/
    d_Z += mfac * d_M + ifac * MCD::Identity(2 * d_N + d_nu, 2 * d_N + d_nu);
    d_slv.compute(d_Z, Eigen::EigenvaluesOnly);
    d_result.row(ctr) = scale * d_slv.eigenvalues().segment(d_N + d_nEig_min, d_nEig_max - d_nEig_min);
    d_det[ctr] = (d_nDet > 0) ? (scale * d_slv.eigenvalues().segment(d_N -(d_nDet / 2), d_nDet)).prod() : 0;
  }

  d_average.resize(2, d_result.cols());
  d_ratios.resize(2, d_result.cols() * (d_result.cols() - 1) / 2);

  bootAver(1000);
  bootRat(1000);
}

void RanMat::bootAver(size_t const nBoot) const
{
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
}

void RanMat::bootRat(size_t const nBoot) const
{
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
}
