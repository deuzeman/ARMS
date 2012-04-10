#include <RanMat.h>

#include <algorithm>
#include <iostream>

RanMat::RanMat(size_t const N, size_t const nu, int const eigMin, int const eigMax)
: d_N(N), d_nu(nu), d_scale(1.0 / static_cast< double >(2 * d_N + d_nu)),
d_eigMin(eigToIndex(eigMin)), d_numEigs(eigMax - eigMin), d_rstream(0),
d_Z(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
d_gamma_5(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
d_A(d_N + d_nu, d_N + d_nu), d_B(d_N, d_N), d_W(d_N + d_nu, d_N),
d_slv(2 * d_N + d_nu), d_result(0)
{
  if ((eigMax - eigMin) < 0)
  {
    std::cerr << "Maximum eigenvalue is smaller than minimum eigenvalue." << std::endl;
    exit(1);
  }
  else if ((N < -2 * eigMin) || N < 2 * eigMax)
  {
    std::cerr << "Requested invalid eigenvalues -- should be equal or less then N/2." << std::endl;
    exit(1);
  }
  
  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);
  
  for (size_t idx = 0; idx < 2 * N + nu; ++idx)
    d_gamma_5(idx, idx) = idx < (N + nu) ? 1.0 : -1.0;
  
  srand(time(0));
  for (size_t ctr = 0; ctr < d_rank; ++ctr)
    rand(); // We fast-forward the stream to be more or less uncorrelated for each process now
  d_rstream.RandomInit(rand()); // Now use the crappy number to seed the Mersenne twister
}

RanMat::~RanMat()
{
  delete[] d_result;
}

void RanMat::calculate(Point const &params, size_t iter)
{
  double const static sqrt1_2 = std::sqrt(2) / 2.0;
  double const static sqrt8   = std::sqrt(8);
  double const static scale2  = d_scale * d_scale;
  double const m  = params.m * scale2;
  double const a6 = params.a6 * scale2 * sqrt8;
  double const a7 = params.a7 * scale2 * sqrt8;
  double const a8 = params.a8 * scale2;
  double const sigma = params.sigma;
  
  iter = (iter / d_nodes) + 1; // Spread the workload 
  iter = ((iter - 1) / 50 + 1) * 50; // Round up to multiples of 50 for jackknife convenience

  delete[] d_result;
  d_result = new double[iter];
  d_samples = iter;
  
  for (size_t ctr = offset; ctr < d_samples; ++ctr)
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
    
    // Compute the eigenvalues and copy the relevant section to the result array
    d_slv.compute(d_Z, Eigen::EigenvaluesOnly);
    for (size_t eig = 0; eig < d_numEigs; ++eig)
      d_result[eig * d_samples + ctr] = ((2 * d_N + d_nu) / sigma) * d_slv.eigenvalues().coeff(d_eigMin + eig);
  }
}

size_t *RanMat::discretize(double const *breaks, int eigMin, size_t const levels, size_t const eigs) const
{
  resDiscrete = new size_t[d_samples * eigs];
  size_t startCol = eigToIndex(eigMin);
  if ((startCol < d_eigMin) || ((startCol + eigs) > (d_eigMin + d_numEigs)))
    exit(1);
  startCol -= d_eigMin;
  for (size_t col = 0; col < eigs; ++col)
    for (size_t row = 0; row < d_samples; ++row)
    {
      double const val = d_result[(startCol + col) * d_samples + row];
      size_t idx = 0;
      while ((idx < levels) && (val > breaks[col * levels + idx]))
        ++idx;
      resDiscrete[col * d_samples + row] = idx;
    }
    return resDiscrete;
}
