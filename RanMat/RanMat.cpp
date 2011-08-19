#include <RanMat.h>

#include <algorithm>
#include <iostream>

RanMat::RanMat(params const *par, int const seed)
  : d_N(par->N),
    d_nu(par->nu),
    d_nEigs(12),
    d_m(par->m),
    d_a6(par->a6),
    d_a7(par->a7),
    d_a8(par->a8),
    d_scale(0.5 / std::sqrt(static_cast< double >(2 * d_N + d_nu))),
    d_Z(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_M(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
    d_A(d_N + d_nu, d_N + d_nu),
    d_B(d_N, d_N),
    d_W(d_N + d_nu, d_N),
    d_slv(2 * d_N + d_nu),
    d_result(2 * d_N + d_nu, d_nEigs),
    d_rstream(seed)

{
  for (int x = 0; x < d_N + d_nu; ++x)
    d_M(x, x) = d_m;
  for (int x = d_N + d_nu; x < (2 * d_N + d_nu); ++x)
    d_M(x, x) = -d_m;
}

void RanMat::changeParams(Eigen::VectorXd const &params)
{
  d_m  = params.coeff(0);
  d_a6 = params.coeff(1);
  d_a7 = params.coeff(2);
  d_a8 = params.coeff(3);
}

void RanMat::calculate(size_t const iter, bool const extend)
{
  size_t offset = 0;

  if (extend)
  {
    offset = d_result.rows();
    d_result.conservativeResize(offset + iter, d_nEigs);
  }
  else
    d_result.resize(iter, d_nEigs);
 
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

    d_Z << d_a8 * d_A, d_W, d_W.adjoint(), d_a8 * d_B;
    
    // The following part now adds the d_W6 and d_W7 terms, as well as the mass term
    for (int x = 0; x < (2 * d_N + d_nu); ++x)
      d_Z(x,x) += (d_rstream.Normal(1.0, d_scale * d_a6)) * d_M(x, x) + d_rstream.Normal(0.0, d_scale * d_a7);

    d_slv.compute(d_Z, Eigen::EigenvaluesOnly);
    d_result.row(ctr) = d_slv.eigenvalues().segment(d_N - (d_nEigs / 2), d_nEigs);
    if ((ctr + 1) % 1000 == 0)
      std::cerr << (ctr + 1) << "\titerations calculated" << std::endl;
  }
}
