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
  d_prec(0.5 * params.prec), 
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

  double result;
  double error = 1.0;
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  
  d_disc.clear();
  while (error > d_prec)
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
    }
    
    result = 0.0;
    for (size_t eig = 0; eig < d_eigs; ++eig)
    {
      double pred = d_disc.average(eig);
      if (d_rank == 0)
      {
	std::cout << "For eigenvalue " << eig << ": measurement = " << d_aver[eig] << ", prediction = " << pred << std::endl;
	std::cout << "                              deviation   =  " << (pred - d_aver[eig]) << std::endl;
      }
      result += std::pow(pred - d_aver[eig], 2.0);
    }
    
    // Now calculate the relative error
    double ave = 0.0;
    double aveSq = 0.0;
    for (size_t bIdx = 0; bIdx < d_blocks; ++bIdx)  
    {
      ave += d_jack[bIdx];
      aveSq += d_jack[bIdx] * d_jack[bIdx];
    }
    ave /= d_blocks;
    aveSq /= d_blocks;
    
    if (d_rank == 0)
      std::cout << aveSq << " - " << ave << "^2 = " << (aveSq - ave * ave) << std::endl;
    error = std::sqrt(aveSq - ave * ave);
    
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << result << " and an error of " << error << '.' << std::endl;

    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((error / d_prec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (error < d_prec && !d_rank)
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
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  double result = 0.0;

  d_disc.clear();
  while (error > d_prec)
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
	double qq = 0.0;
        for (size_t samp = 0; samp < d_levels; ++samp)
        {
          double cc = std::abs(d_disc(eig, samp, blockIdx) - samp * d_inc);
          qq = std::max(qq, cc);
        }
        d_jack[blockIdx] += qq;
      }
      d_jack[blockIdx] /= d_eigs;
    }
    
    result = 0.0;
    for (size_t eig = 0; eig < d_eigs; ++eig)
    {
      double qq = 0.0;
      for (size_t samp = 1; samp < d_levels; ++samp)
      {
        double cc = std::abs(d_disc(eig, samp) - samp * d_inc);
        qq = std::max(qq, cc);
      }
      result += qq;
    }
    result /= d_eigs;
      
    // Now calculate the relative error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    double ave = 0.0;
    double aveSq = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
    {
      ave += d_jack[idx];
      aveSq += d_jack[idx] + d_jack[idx];
    }
    ave /= d_blocks;
    aveSq /= d_blocks;
    error = std::sqrt(aveSq - ave * ave);
    
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << result << " and an error of " << error << '.' << std::endl;
    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((error / d_prec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (error < d_prec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  
  return result;
}

Comparator::~Comparator()
{
  delete[] d_jack;
}
