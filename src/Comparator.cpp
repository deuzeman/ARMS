#include <algorithm>
#include <numeric>
#include <Comparator.h>
#include <Log.h>

Comparator::Comparator(Data &data, Params &params, int type)
: d_ranmat(params.N, params.nu, data.minEv(), data.maxEv()), 
  d_aver(data.average()),
  d_av_err(data.error()),
  d_levels(data.numSamples() + 1), 
  d_inc(1.0 / data.numSamples()), 
  d_eigs(data.numCols()), 
  d_minEv(data.minEv()),
  d_blocks(params.blocks),
  d_prec_a(0.5 * params.prec_a),
  d_prec_k(0.5 * params.prec_k),
  d_type(type),
  d_disc(data.flatPerColumn(), data.numSamples(), d_eigs, d_blocks),
  d_jack(new double[d_blocks]),
  d_cumdist(new double[d_levels])
{
  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);
}

double Comparator::deviation(Point const &point)
{
  if (Log::ionode)
    log() << "  >>  Calculating deviation for " << point << std::endl;

  if (Log::ionode)
  {
    if (d_type == AVE)
      log() << "  >>  Calculating peak average chi squared." << std::endl;
    else
      log() << "  >>  Calculating Kolmogorov-Smirnov D." << std::endl;
  }
  
  // Set up some data structure for the algorithms
  double error = 1.0;
  size_t needed = roundToBlocks(2 * d_levels); // We don't want the resolution here to be an issue
  
  size_t samples = 0;
  double result = 0.0;

  d_disc.clear();

  double prec = (d_type == AVE) ? d_prec_a : d_prec_k;
  while ((error > prec) && (samples < 1000000))
  {
    if (Log::ionode)
      log() << "  >>  Requesting " << needed << " samples." << std::endl;

    // Calculate as much data as we think we will require
    samples += d_nodes * d_ranmat.calculate(point, needed);
    d_disc.calculate(d_ranmat);

    // Calculate the main result
    result = 0.0;
    if (d_type == AVE)
      for (size_t eig = 0; eig < d_eigs; ++eig)
      {
        double chi = (d_disc.average(eig) - d_aver[eig]) / d_av_err[eig];
        result += chi * chi;
      }
    else
      for (size_t eig = 0; eig < d_eigs; ++eig)
        for (size_t samp = 0; samp < d_levels; ++samp)
        {
          result = std::max(result, std::abs(d_disc(eig, samp) - samp * d_inc));
//          std::cerr << "[DEBUG] |" << d_disc(eig, samp) << " - " << (samp * d_inc) << "| = " << std::abs(d_disc(eig, samp) - samp * d_inc) << std::endl;
	}

    // Extract the jackknife samples
    for (size_t block = 0; block < d_blocks; ++block)
    {
      d_jack[block] = 0.0;
      for (size_t eig = 0; eig < d_eigs; ++eig)
      {
        if (d_type == AVE)
        {
          double chi = (d_disc.average(eig, block) - d_aver[eig]) / d_av_err[eig];
          d_jack[block] += chi * chi;
        }
        else
        {
          for (size_t samp = 0; samp < d_levels; ++samp)
            d_jack[block] = std::max(d_jack[block], std::abs(d_disc(eig, samp, block) - samp * d_inc));
        }
      }
    }
          
    // Now calculate the absolute error
    // Remember: this is a jackknife, so we *sum* over differences squared and do the rescaling
    error = 0.0;
    for (size_t idx = 0; idx < d_blocks; ++idx)
    {
      error += (d_jack[idx] - result) * (d_jack[idx] - result);
    }
    error *= (d_blocks - 1.0) / d_blocks;
    error = std::sqrt(error);

    if (result > 1.0) // In these cases, we want to make the error relative
      error /= result;

    if (Log::ionode)
      log() << "  >>  With a total of " << samples << " samples, obtained a value of " << result << " and an error of " << error << '.' << std::endl;
    
    // We'll add a minimum and maximum number of iterations
    // To avoid waiting forever for the 1M measurements
    // or watching the thing skip around in increments of 50.
    needed = std::min(std::max(roundToBlocks(static_cast< size_t >(std::pow((error / prec), 2.0) * samples)), 
                                             static_cast< size_t >(1000)), static_cast< size_t >(50000));
    if (error < prec && !d_rank)
      log() << "  >>  This is sufficient for the currently needed precision.\n" << std::endl;
  }
  if (samples > 1000000)
    if (Log::ionode)
      log() << " >> [WARNING] Sample number exceeding 1M!\n" << std::endl;
  return result;
}

Comparator::~Comparator()
{
  delete[] d_jack;
  delete[] d_aver;
  delete[] d_av_err;
}

