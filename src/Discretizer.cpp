#include <Discretizer.h>
#include <Log.h>
#include <numeric>

Discretizer::Discretizer(double const *breaks, unsigned long numBreaks, unsigned long numEigs, unsigned long blocks)
 : d_numEigs(numEigs), d_numBlocks(blocks), d_numLevels(numBreaks + 1), d_sampTotal(0), d_breaks(breaks)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);

   d_sum_blocks = new double[d_numEigs * d_numBlocks];
   d_mean_blocks = new double[d_numEigs * d_numBlocks];

   d_hist_blocks = new unsigned long[d_numEigs * d_numLevels * d_numBlocks];
   d_cum_blocks = new double[d_numEigs * d_numLevels * d_numBlocks];

   clear();
}

void Discretizer::clear()
{
  d_sampTotal = 0;
  std::fill_n(d_sum_blocks, d_numEigs * d_numBlocks, 0.0);
  std::fill_n(d_hist_blocks, d_numEigs * d_numLevels * d_numBlocks, 0.0);
}

void Discretizer::calculate(RanMat const &ranmat)
{
  double const *res = ranmat.result();

  // Sum the data block-by-block
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
    for (unsigned long samp = 0; samp < ranmat.numSamples(); ++samp)
      d_sum_blocks[(samp % d_numBlocks) * d_numEigs + eig] += res[eig * ranmat.numSamples() + samp];
  // This is all that is needed for the calculation of averages

  // Gather the data as a histogram
  unsigned long *data = new unsigned long[ranmat.numSamples() * d_numEigs];
  std::fill_n(data, ranmat.numSamples() * d_numEigs, 0);

  // Loop over the data and get the binning numbers
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
  {
    unsigned long *pd = data + eig * ranmat.numSamples();
    double const *pb = d_breaks + eig * (d_numLevels - 1);
    double const *pr = res + eig * ranmat.numSamples();
    for (unsigned long samp = 0; samp < ranmat.numSamples(); ++samp)
      while ((pd[samp] < (d_numLevels - 1)) && pr[samp] > pb[pd[samp]])
        ++pd[samp];
  }

  // Now go over the binning numbers and produce bin counts
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
  {
    unsigned long *pd = data + eig * ranmat.numSamples();
    for (unsigned long samp = 0; samp < ranmat.numSamples(); ++samp)
      ++d_hist_blocks[((samp % d_numBlocks)  * d_numEigs + eig) * d_numLevels + pd[samp]];
  }
  delete[] data;

  d_sampTotal += d_nodes * ranmat.numSamples();

  // We have now gathered the overall data per block per node.
  // These are just sums, so we don't need to do any rescaling.
  // Now we want to produce averages globally.

  // Aggregate the sums globally and calculate the averages per block
  double avfac = static_cast< double >(d_numBlocks) / d_sampTotal;
  MPI_Allreduce(d_sum_blocks, d_mean_blocks, d_numEigs * d_numBlocks, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  for (unsigned long idx = 0; idx < d_numEigs * d_numBlocks; ++idx)
    d_mean_blocks[idx] *= avfac;

  // Aggregate the histograms globally
  unsigned long *glob_hist_blocks = new unsigned long[d_numEigs * d_numLevels * d_numBlocks];
  MPI_Allreduce(d_hist_blocks, glob_hist_blocks, d_numEigs * d_numLevels * d_numBlocks, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned long idx = 0; idx < d_numEigs * d_numLevels * d_numBlocks; ++idx)
    d_cum_blocks[idx] = avfac * glob_hist_blocks[idx];
  delete[] glob_hist_blocks;

  // Now create a cumulative distribution from the global histograms
  for (unsigned long block = 0; block < d_numBlocks; ++block)
    for (unsigned long eig = 0; eig < d_numEigs; ++eig)
    {
      for (unsigned long level = 1; level < d_numLevels; ++level)
        d_cum_blocks[(block * d_numEigs + eig) * d_numLevels + level] += d_cum_blocks[(block * d_numEigs + eig) * d_numLevels + level - 1];
    }
}

double Discretizer::operator()(unsigned long const &eig, unsigned long const &level) const
{
  double res = 0.0;
  for (unsigned long bIdx = 0; bIdx < d_numBlocks; ++bIdx)
    res += d_cum_blocks[(bIdx * d_numEigs + eig) * d_numLevels + level];
  return (res / d_numBlocks);
}

double Discretizer::average(unsigned long const &eig) const
{
  double res = 0.0;
  for (unsigned long bIdx = 0; bIdx < d_numBlocks; ++bIdx)
    res += d_mean_blocks[bIdx * d_numEigs + eig];
  return (res / d_numBlocks);
}


Discretizer::~Discretizer()
{
  delete[] d_sum_blocks;
  delete[] d_mean_blocks;
  
  delete[] d_hist_blocks;
  delete[] d_cum_blocks;
}
