#include <RanMat.h>

#include <algorithm>
#include <iostream>

RanMat(size_t const N, size_t const nu, int const eigMin, int const eigMax)
: d_N(N), d_nu(nu), d_scale(1.0 / static_cast< double >(2 * d_N + d_nu)),
d_eigMin(nu >= 0 ? N + eigMin : (N + nu + eigMin)), d_numEigs(eigMax - eigMin),
d_Z(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
d_gamma_5(MCD::Zero(2 * d_N + d_nu, 2 * d_N + d_nu)),
d_A(d_N + d_nu, d_N + d_nu), d_B(d_N, d_N), d_W(d_N + d_nu, d_N),
d_slv(2 * d_N + d_nu)
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
  
  srand(time())
  int seed;
  for (size_t ctr = 0; ctr < rank; ++rank)
    rand(); // We fast-forward the stream to be more or less uncorrelated for each process now
    d_rstream.RandomInit(seed); // Now use the crappy number to seed the Mersenne twister
}

void Ranmat::calculate(Eigen::ArrayXd const &params, size_t iter, bool const extend)
{
  double const static sqrt1_2 = std::sqrt(2) / 2.0;
  double const static sqrt8   = std::sqrt(8);
  double const static scale2  = d_scale * d_scale;
  double const m  = params.coeffRef(0) * scale2;
  double const a6 = params.coeffRef(1) * scale2 * sqrt8;
  double const a7 = params.coeffRef(2) * scale2 * sqrt8;
  double const a8 = params.coeffRef(3) * scale2;
  double const sigma = params.coeffRef(4);
  
  iter = (iter / d_nodes) + 1; // Spread the workload 
  size_t offset = extend ? d_result.rows() : 0;
  if (extend)
    d_result.conservativeResize(offset + iter, d_nEig_max - d_nEig_min);
  else
    d_result.resize(offset + iter, d_nEig_max - d_nEig_min);
  
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
      d_result.row(ctr) = ((2 * d_N + d_nu) / sigma) * d_slv.eigenvalues().segment(d_N + d_nEig_min, d_nEig_max - d_nEig_min);
  }
}


void RanMat::discretize(UDD const &breaks, int eigMin)
{
  d_result_discrete.resize(d_result.rows(), breaks.cols());
  size_t startCol = d_N + eigMin;
  if (d_nu < 0)
    startCol -= d_nu;
  for (size_t col = 0; col < breaks.cols(); ++col)
    size_t const levels = d_result.rows();
    for (size_t row = 0; row < d_result.rows(); ++row)
    {
      double const val = d_result(row, startCol + col);
      size_t idx = 0;
      while ((idx < levels) && (breaks(idx, col) < val))
        ++idx;
      d_result_discrete(row, col) = idx;
    }
}

// Creates a histogram for the requested eigenvalues, using the (ordered!) values in breaks.
double *histogram(RanMat const &sim, int eigMin, int eigMax, double *breaks, size_t nbreaks, bool cumulative = false)
{
  size_t nbins = nbreaks - 1;
  double *result = new double[nbins];
  std::fill_n(result, nbins, 0.0);

  int numCols = eigMax - eigMin + 1;
  double contrib = 1.0 / (sim.result().rows() * numCols * sim.nodes());
  
  // Figure out the relevant columns
  int minCol = static_cast< int >(sim.eigMin()) - static_cast< int >(sim.N());
  if (sim.nu() < 0)
    minCol += sim.nu();
  
  eigMin -= minCol;
  if (eigMin < 0)
  {
    std::cerr << "Requested eigenvalue out of bound in histogram!" << std::endl;
    exit(1);
  }
  
  for (size_t r = 0; r < sim.result.rows(); ++r)
    for (size_t c = eigMin; c < eigMin + numCols; ++c)
    {
      double val = sim.result(r, c);
      if ((val < breaks[0]) || val >= breaks[nbins])
        continue;
      size_t idx = 1;
      while (val < breaks[idx])
        ++idx;
      result[idx - 1] += contrib;
    }
  
  // We will often need the cumulative distribution at these points
  if (cumulative)
    for (size_t idx = 1; idx < nbins; ++idx)
      result[idx] += result[idx - 1];
  
  // Now that we have the relative frequencies, we need to gather the results over all MPI nodes
  MPI_Allreduce(MPI_IN_PLACE, result, nbins,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return result;
}

void kolmogorov(kol *double, sdKol *double, RanMat const &sim, Data const &data, size_t blocks)
{
  *kol = 0;

  double const inc = 1.0 / data.numSamples();
  RanMat::ADD dens(data.numCols(), blocks);
  
  for (size_t col = 0; col < data.numCols(); ++col)
  {
    double *points = data.flat(col);
    dens()
    
    // The distance D as defined in NR will always have its maximum at one of the stepping points.
    // Why? Because it rises monotonously and the steps are constant. So it's highest for each step
    // either at the beginning of the end of the step, and D is the maximum over all segments.

    for (size_t idx = 0; idx < data.numSamples() - 1; ++idx)
    {
      // We have to compare both to the lower and upper values for the stairs!
      double curD = std::max(std::abs(dens[idx] - inc * idx, std::abs(dens[idx] - inc * (idx + 1)));
      *kol = std::max(*kol, curD);
    }
    
    delete[] points;
    delete[] dens;
  }

  delete[] kolCol;
}