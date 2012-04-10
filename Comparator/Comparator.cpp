#include <algorithm.h>
#include <numeric.h>
#include <Comparator.h>

double Comparator::kolmogorov(Point const &point)
{ 
  // Set up some data structure for the algorithms
  double error = 1.0;
  double const rescale = static_cast< double >(d_blocks) / (d_blocks - 1.0);
  size_t needed = roundToBlocks(1000); // Sounds like a decent number to start with?
  
  double *dists = new double*[d_blocks * (d_levels + 1) * d_eigs];
  double *full  = new double[(d_levels + 1) * d_eigs];
  double *subt  = new double[(d_levels + 1) * d_eigs];
  double *jack  = new double[d_blocks];
  double *cum   = new double[(d_levels + 1) * d_eigs];
  
  std::fill_n(dists, d_blocks * (d_levels + 1) * d_eigs, 0.0);
  std::fill_n(full, (d_levels + 1) * d_eigs);
  
  size_t samples = 0;
  double result = 0.0;
  
  while (error > d_relprec)
  {
    // Calculate as much data as we think we will require
    d_ranmat.calculate(point, needed);
    samples += needed;
    double contrib = 1.0 / samples;
    
    // Discretize the newly calculated data and make a histogram from it.
    // Note that we loop over the numbers of blocks so we can take a jackknife
    // error estimate later on.
    size_t *dres = d_ranmat.discretize(breaks, d_minEv, d_levels, d_eigs);
    size_t blckCtr = 0;
    for (size_t col = 0; col < d_eigs; ++col)
      for (size_t row = 0; row < (needed / d_nodes); ++row)
      {
        full[col * d_levels + dres[col * (needed / d_nodes) + row]] += contrib;
        dists[(blckCtr * (d_levels + 1) * d_eigs) +  col * d_levels + dres[col * (needed / d_nodes) + row]] += contrib;
        ++blckCtr;
        blckCtr %= d_blocks;
      }
    delete[] dres;
    
    // Aggregate the data to make cumulative a distribution -- we will do this below for the blocks
    std::partial_sum(full, full + (d_levels + 1) * d_eigs, cum);
    std::copy(cum, (d_levels + 1) * d_eigs, full);
    
    // Perform a jackknife error estimate -- calculate all the samples
    for (size_t blockIdx = 0; blockIdx < d_blocks; ++blockIdx)
    {
      // Calculate the jackknife sample -- first taking the cumulative of the separate blocks
      // We do that here, because there's no point in copying those values back.
      std::partial_sum(dists + blockIdx * (d_levels + 1) * d_eigs, dists + (blockIdx + 1) * (d_levels + 1) * d_eigs, cum);
      for (size_t idx = 0; idx < (d_levels + 1) * d_eigs; ++idx)
        subt[idx] = rescale * (full[idx] - cum[idx]);

      // Calculate the cumulative sum
      MPI_Allreduce(MPI_IN_PLACE, subt, (d_levels + 1) * d_eigs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
      // Calculate the maximum deviation
      // NOTE Check -- we probably need to take into account both sides of the step function!
      jack[blockIdx] = 0.0;
      for (size_t col = 0; col < d_eigs; ++col)
        for (size_t row = 1; row < d_levels; ++row)
          jack[blockIdx] = std::max(jack[blockIdx], std::abs(cum[col * (d_levels + 1) + row] - col * inc));
    }
    
    // At this point, jack contains an array of jackknife averages. Now we need the proper estimate.
    MPI_Allreduce(MPI_IN_PLACE, full, (d_levels + 1) * d_eigs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    result = 0.0;
    for (size_t col = 0; col < d_eigs; ++col)
      for (size_t row = 1; row < d_levels; ++row)
        result = std::max(result, std::abs(full[col * (d_levels + 1) + row] - col * inc));
      
    // Now calculate the relative error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    double ave = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
      ave += jack[idx];
    ave /= d_blocks;
    
    error = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
      error += (jack[idx] - ave) * (jack[idx] - ave);
    error *= rescale;
    error /= result;
  }
  
  delete[] dists;
  delete[] full;
  delete[] subt;
  delete[] jack;
  delete[] cum;
  
  return result;
}
