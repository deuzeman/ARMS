#include <Discretizer.h>

#include <numeric>

Discretizer::Discretizer(double *breaks, unsigned long numBreaks, unsigned long numEigs, unsigned long blocks)
 : d_numEigs(numEigs), d_numBlocks(blocks), d_numLevels(numBreaks + 1), d_sampTotal(0), d_breaks(breaks)
{
   MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);
   
   d_sum_blocks = new double[d_numEigs * d_numBlocks];
   d_mean_blocks = new double[d_numEigs * d_numBlocks];
   
   d_hist_blocks = new unsigned long[d_numLevels * d_numEigs * d_numBlocks];
   d_cum_blocks = new double[d_numLevels * d_numEigs * d_numBlocks];

   d_mean_total = new double[d_numEigs];
   d_cum_total = new double[d_numLevels * d_numEigs];

   clear();
}

void Discretizer::clear()
{
  d_sampTotal = 0;
  std::fill_n(d_sum_blocks, d_numEigs * d_numBlocks, 0.0);
  std::fill_n(d_hist_blocks, d_numLevels * d_numEigs * d_numBlocks, 0.0);
}
 
void Discretizer::calculate(RanMat const &ranmat)
{
  double const *res = ranmat.result();
  unsigned long rms = ranmat.numSamples();

  // Sum the data block-by-block
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
    for (unsigned long samp = 0; samp < rms; ++samp)
      d_sum_blocks[(samp % d_numBlocks) * d_numEigs + eig] += res[eig * rms + samp];
  // This is all that is needed for the calculation of averages
  
  // Gather the data as a histogram
  unsigned long *data = new unsigned long[rms * d_numEigs];
  std::fill_n(data, rms * d_numEigs, 0);

  // Loop over the data and get the bin indices for all RM samples
  // Should now work for subsequent arrays of sorted values per eigenvalue
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
  {
    unsigned long *data_col = data + eig * rms;
    double const *res_col = res + eig * rms;
    for (unsigned long samp = 0; samp < rms; ++samp)
      while ((data_col[samp] < d_numLevels) && res_col[samp] > d_breaks[(d_numLevels - 1) + data_col[samp]])
        ++data_col[samp];
  }

  // Now go over the binning numbers and produce bin counts
  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
  {
    unsigned long *data_col = data + eig * rms;
    // Note that we treat the jackknife blocks as the outermost index here
    for (unsigned long samp = 0; samp < rms; ++samp)
      ++d_hist_blocks[(samp % d_numBlocks) * d_numLevels * d_numEigs + d_numLevels * eig + data_col[samp]];
  }

  // The binning data can be removed, it was just intermediate notation
  delete[] data;
  
  // We have now gathered the overall data per block per node.
  // These are just sums, so we don't need to do any rescaling.
  // Now we want to produce averages globally.
  
  // Aggregate the sums globally  
  MPI_Allreduce(d_sum_blocks, d_mean_blocks, d_numEigs * d_numBlocks, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Calculate the averages per block
  d_sampTotal += d_nodes * ranmat.numSamples();
  double avfac = static_cast< double >(d_numBlocks) / d_sampTotal;

  for (unsigned long eig = 0; eig < d_numEigs; ++eig)
  {
    d_mean_total[eig] = 0;
    for (unsigned long block = 0; block < d_numBlocks; ++block)
    {
      d_mean_blocks[block * d_numEigs + eig] *= avfac;
      d_mean_total[eig] += d_mean_blocks[block * d_numEigs + eig] / d_numBlocks;
    }
  }
  
  // Aggregate the histograms globally
  unsigned long *glob_hist_blocks = new unsigned long[d_numLevels * d_numEigs * d_numBlocks];
  MPI_Allreduce(d_hist_blocks, glob_hist_blocks, d_numLevels * d_numEigs * d_numBlocks, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  double cumfac = static_cast< double >(d_numBlocks) / d_sampTotal;
  for (unsigned long idx = 0; idx < d_numLevels * d_numEigs * d_numBlocks; ++idx)
    d_cum_blocks[idx] = cumfac * glob_hist_blocks[idx];
  delete[] glob_hist_blocks;
  
  // Now create a cumulative distribution from the global histograms
  for (unsigned long block = 0; block < d_numBlocks; ++block)
    for (unsigned long eig = 0; eig < d_numEigs; ++eig)
      for (unsigned long level = 1; level < d_numLevels; ++level)
        d_cum_blocks[block * d_numLevels * d_numEigs + d_numLevels * eig + level] += d_cum_blocks[block * d_numLevels * d_numEigs + d_numLevels * eig + level - 1];

  // Produce one sum distribution from the separate jackknife blocks
  std::fill_n(d_cum_total, d_numLevels * d_numEigs, 0.0);
  for (unsigned long block = 0; block < d_numBlocks; ++block)
    for (unsigned long idx = 0; idx < d_numLevels * d_numEigs; ++idx)
      d_cum_total[idx] += d_cum_blocks[block * d_numLevels * d_numEigs + idx];
}

Discretizer::~Discretizer()
{
  delete[] d_breaks;
  
  delete[] d_sum_blocks;
  delete[] d_mean_blocks;
  
  delete[] d_hist_blocks;
  delete[] d_cum_blocks;
}
