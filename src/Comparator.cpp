#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
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
  d_point = &point;

  if (Log::ionode)
    log() << "  >>  Calculating eigenvalue average deviation for " << point << std::endl;

  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  d_error = 1.0;  
  d_disc.clear();
  while ((d_error > d_prec) && samples < 2500000)
  {
    if (Log::ionode)
      log() << "  >>  Requesting " << needed << " samples." << std::endl;
      // Calculate as much data as we think we will require

    d_ranmat.calculate(point, needed);
    samples += needed;
    d_disc.calculate(d_ranmat);
    
    for (size_t block = 0; block < d_blocks; ++block)
    {
      d_jack[block] = 0.0;
      for (size_t eig = 0; eig < d_eigs; ++eig)
      {
        double pred = d_disc.average(eig, block);
        d_jack[block] += std::pow((pred - d_aver[eig]) / d_aver[eig], 2.0);
      }
    }
    
    d_result = 0.0;
    for (size_t eig = 0; eig < d_eigs; ++eig)
    {
      double pred = d_disc.average(eig);
      d_result += std::pow((pred - d_aver[eig]) / d_aver[eig], 2.0);
    }
    
    // Now calculate the relative d_error
    double ave = 0.0;
    double aveSq = 0.0;
    for (size_t block = 0; block < d_blocks; ++block)  
    {
      ave += d_jack[block];
      aveSq += d_jack[block] * d_jack[block];
    }
    ave /= d_blocks;
    aveSq /= d_blocks;
    
    d_error = std::sqrt((aveSq - ave * ave) / d_blocks);
    
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << d_result << " and an d_error of " << d_error << '.' << std::endl;

    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((d_error / d_prec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (d_error < d_prec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  if (samples > 1000000)
    if (Log::ionode)
      log() << " >> [WARNING] Sample number exceeding 1M!\n" << std::endl;
  return d_result;
}

double Comparator::kolmogorov(Point const &point)
{
  d_point = &point;
  
  if (Log::ionode)
    log() << "  >>  Calculating Kolmogorov-Smirnov D for " << point << std::endl;
  
  // Set up some data structure for the algorithms
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  d_error = 1.0;

  d_disc.clear();
  while (d_error > d_prec && samples < 1000000)
  {
    if (Log::ionode)
      log() << "  >>  Requesting " << needed << " samples." << std::endl;
    // Calculate as much data as we think we will require
    d_ranmat.calculate(point, needed);
    samples += needed;
    d_disc.calculate(d_ranmat);

    for (size_t block = 0; block < d_blocks; ++block)
    {
      d_jack[block] = 0.0;
      for (size_t eig = 0; eig < d_eigs; ++eig)
        for (size_t samp = 0; samp < d_levels; ++samp)
          d_jack[block] = std::max(d_jack[block], std::abs(d_disc(eig, samp, block) - samp * d_inc));
    }
    
    d_result = 0.0;
    for (size_t eig = 0; eig < d_eigs; ++eig)
      for (size_t samp = 1; samp < d_levels; ++samp)
        d_result = std::max(d_result, std::abs(d_disc(eig, samp) - samp * d_inc));
      
    // Now calculate the relative error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    double ave = 0.0;
    double aveSq = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
    {
      ave += d_jack[idx];
      aveSq += d_jack[idx] * d_jack[idx];
    }
    ave /= d_blocks;
    aveSq /= d_blocks;
    d_error = std::sqrt((aveSq - ave * ave) / d_blocks);
    
    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << d_result << " and an error of " << d_error << '.' << std::endl;
    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = roundToBlocks(std::min(std::max(static_cast< size_t >(std::pow((d_error / d_prec), 2.0) * samples), 1000ul), 50000ul));
    if (d_error < d_prec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  if (samples > 1000000)
    if (Log::ionode)
      log() << " >> [WARNING] Sample number exceeding 1M!\n" << std::endl;
  return d_result;
}

Comparator::~Comparator()
{
  delete[] d_jack;
}

void Comparator::writeCumulative() const
{
  if (!Log::ionode)
    return;

  std::cout << "# Cumulative distribution for " << *d_point << ",\n"
            << "# with a value of " << d_result << " +/- " << d_error << std::endl;
  std::cout << "Eigenvalue" << std::setw(15) << std::right << "break" << std::setw(15) << std::right << "density" << std::endl;
  for (int eig = 0; eig < d_eigs; ++eig)
  {
    int hEig = d_minEv + eig;
    hEig = hEig < 0 ? hEig : hEig + 1;
    for (unsigned int level = 1; level < d_levels; ++level)
      std::cout << std::setw(10) << std::right << hEig 
                << std::setw(15) << std::right << d_breaks[eig * (d_levels - 1) + level] 
                << std::setw(15) << std::right << d_disc(eig, level) << std::endl;
  }
}
