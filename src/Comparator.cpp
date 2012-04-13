#include <algorithm>
#include <numeric>
#include <Comparator.h>
#include <Log.h>

Comparator::Comparator(Data &data, Params &params)
: d_ranmat(params.N, params.nu, data.minEv(), data.maxEv()), 
d_breaks(data.flatPerColumn()),
d_aver(data.average()),
d_levels(data.numSamples() + 1), 
d_inc(1.0 / data.numSamples()), 
d_eigs(data.numCols()), 
d_minEv(data.minEv()),
d_blocks(params.blocks),
d_relprec(0.1 * params.prec), 
d_disc(d_breaks, data.numSamples(), d_eigs, d_blocks),
d_jack(new double[d_blocks])
{
  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);
}

double Comparator::averages(Point const &point)
{
  if (Log::ionode)
    log() << "  >>  Calculating eigenvalue average deviation for " << point << std::endl;
  
  double result = 0.0;
  double error = 1.0;
  double const rescale = static_cast< double >(d_blocks) / (d_blocks - 1.0);
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  
  d_disc.clear();
  while (error > d_relprec)
  {
    if (Log::ionode)
      log() << "  >>  Requesting " << needed << " samples." << std::endl;
      // Calculate as much data as we think we will require

    d_ranmat.calculate(point, needed);
    samples += needed;
    d_disc.calculate(d_ranmat);
    
    for (size_t blockIdx = 0; blockIdx < d_blocks; ++blockIdx)
    {
      d_jack[blockIdx] = 0.0;
      for (size_t eig = 0; eig < d_eigs; ++eig)
      {
        double pred = d_disc.average(eig, blockIdx);
        d_jack[blockIdx] += std::pow(pred - d_aver[eig], 2.0);
      }
      std::cerr << "Block " << blockIdx << ":  " << d_jack[blockIdx] << std::endl;
    }
    for (size_t eig = 0; eig < d_eigs; ++eig)
    {
      double pred = d_disc.average(eig);
      result += std::pow(pred - d_aver[eig], 2.0);
    }
    std::cerr << "Result: " << result << std::endl;
    
    // Now calculate the relative error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    double ave = 0.0;
    for (size_t bIdx = 0; bIdx < d_blocks; ++bIdx)  
      ave += d_jack[bIdx];
    ave /= d_blocks;
    std::cerr << "Average: " << ave << std::endl;
    
    error = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
      error += (d_jack[idx] - ave) * (d_jack[idx] - ave);
    error = std::sqrt(error * rescale);
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << result << " and an error of " << error << '.' << std::endl;
    error /= result;
    if (Log::ionode)
      log() << "  >>  That implies a relative error of " << error << '.' << std::endl;

    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((error / d_relprec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (error < d_relprec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  return result;
}

double Comparator::kolmogorov(Point const &point)
{
  if (Log::ionode)
    log() << "  >>  Calculating Kolmogorov-Smirnov D for " << point << std::endl;
  
  // Set up some data structure for the algorithms
  double error = 1.0;
  double const rescale = static_cast< double >(d_blocks) / (d_blocks - 1.0);
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  double result = 0.0;

  d_disc.clear();
  while (error > d_relprec)
  {
    if (Log::ionode)
      log() << "  >>  Requesting " << needed << " samples." << std::endl;
    // Calculate as much data as we think we will require
    d_ranmat.calculate(point, needed);
    samples += needed;
    d_disc.calculate(d_ranmat);

    for (size_t blockIdx = 0; blockIdx < d_blocks; ++blockIdx)
    {
      d_jack[blockIdx] = 0.0;
      for (size_t eig = 0; eig < d_eigs; ++eig)
        for (size_t samp = 0; samp < d_levels; ++samp)
        {
          double cc = std::abs(d_disc(eig, samp, blockIdx) - samp * d_inc);
            d_jack[blockIdx] += cc;
        }
      d_jack[blockIdx] /= (d_levels * d_eigs);
    }
    
    result = 0.0;
    for (size_t eig = 0; eig < d_eigs; ++eig)
      for (size_t samp = 1; samp < d_levels; ++samp)
      {
        double cc = std::abs(d_disc(eig, samp) - samp * d_inc);
          result += cc;
      }
      result /= (d_levels * d_eigs);
      
    // Now calculate the relative error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    double ave = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
      ave += d_jack[idx];
    ave /= d_blocks;
    
    error = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
      error += (d_jack[idx] - ave) * (d_jack[idx] - ave);
    error = std::sqrt(error * rescale);
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << result << " and an error of " << error << '.' << std::endl;
    error /= result;
    if (Log::ionode)
      log() << "  >>  That implies a relative error of " << error << '.' << std::endl;

    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((error / d_relprec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (error < d_relprec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  
  return result;
}

Comparator::~Comparator()
{
  delete[] d_jack;
}
