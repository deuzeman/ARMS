#include <Discretizer.h>

#include <numeric>

Discretizer::Discretizer(double const *breaks, size_t numBreaks, size_t numEigs, size_t blocks)
 : d_numEigs(numEigs), d_levels(numBreaks + 1), d_numBlocks(blocks), d_numBreaks(numBreaks), d_breaks(breaks), d_data(0), d_blockFac(d_numBlocks / (d_numBlocks - 1.0))
{
   MPI_Comm_size(MPI_COMM_WORLD, &d_nodes);
   d_acc = new double[d_numBlocks];
   d_hist = new double[d_levels * d_numBlocks];
   d_cum = new double[d_numEigs * d_levels];
   d_block = new double[d_numEigs * d_levels * d_numBlocks];
   d_avBlocks = new double[d_numEigs * d_numBlocks];
   
   clear();
}

void Discretizer::clear()
{
  std::fill_n(d_cum, d_numEigs * d_levels, 0.0);
  std::fill_n(d_block, d_numEigs * d_levels * d_numBlocks, 0.0);
  std::fill_n(d_avBlocks, d_numEigs * d_numBlocks, 0.0);
  
  d_sampTotal = 0;
}
 
void Discretizer::calculate(RanMat const &ranmat)
{
  d_levels = d_numBreaks + 1; 
  
  double const *res = ranmat.result();
  d_data = new size_t[ranmat.numSamples() * d_numEigs];
  std::fill_n(d_data, ranmat.numSamples() * d_numEigs, 0);
  
  for (size_t eig = 0; eig < d_numEigs; ++eig)
  {
    size_t *pd = d_data + eig * ranmat.numSamples();
    double const *pb = d_breaks + eig * d_numBreaks;
    for (size_t samp = 0; samp < ranmat.numSamples(); ++samp)
      while ((pd[samp] < d_numBreaks) && res[samp] > pb[pd[samp]])
        ++pd[samp];
  }

  // Check if we are extending our data here.
  // If so, fix the normalization of the data here.
  if (d_sampTotal > 0)
  {
    double refac = static_cast< double >(d_sampTotal) / (d_sampTotal + d_nodes * ranmat.numSamples());
    // d_cum will be overwritten later on, so we don't need to bother rescaling
    std::fill_n(d_cum, d_levels * d_numEigs, 0.0);
    for (size_t idx = 0; idx < d_numEigs * d_levels * d_numBlocks; ++idx)
      d_block[idx] *= refac;
    for (size_t idx = 0; idx < d_numEigs * d_numBlocks; ++idx)
      d_avBlocks[idx] *= refac;
  }
  d_sampTotal += d_nodes * ranmat.numSamples();
  
  double const inc = 1.0 / d_sampTotal;
  
  size_t const bStep = d_numEigs * d_levels;

  size_t const bSize = ranmat.numSamples() / d_numBlocks;
  double const avfac = static_cast< double >(ranmat.numSamples()) / (d_sampTotal * bSize);
  
  for (size_t eig = 0; eig < d_numEigs; ++eig)
  {
    std::fill_n(d_hist, d_levels * d_numBlocks, 0.0);

    size_t *pd = d_data + eig * ranmat.numSamples();
    double *pc = d_cum + eig * d_levels;
    
    for (size_t samp = 0; samp < ranmat.numSamples(); ++samp)
    {
      d_hist[(samp % d_numBlocks) * d_levels + pd[samp]] += inc;
    }
    
    // If we're doing an MPI run, this is the point to aggregate those results
    MPI_Allreduce(MPI_IN_PLACE, d_hist, d_levels * d_numBlocks, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // In the remainder, all data are already global sums, so we don't need to think about th5250027is further.
    std::fill_n(d_acc, d_numBlocks, 0.0);
    for (size_t lev = 0; lev < d_levels; ++lev)
    {
      for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
      {
        d_acc[bIdx] += d_hist[bIdx * d_levels + lev];
        d_block[(bIdx * d_numEigs + eig) * d_levels + lev] += d_acc[bIdx];
        pc[lev] += d_block[(bIdx * d_numEigs + eig) * d_levels + lev];
      }
    }
    
    for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
    {
      d_avBlocks[eig * d_numBlocks + bIdx] += 
                avfac * std::accumulate(res + eig * ranmat.numSamples() + bIdx * bSize, res + eig * ranmat.numSamples() + (bIdx + 1) * bSize, 0.0);
    }
    
    // If we're doing an MPI run, aggregate the results on the averages...
  }
  MPI_Allreduce(MPI_IN_PLACE, d_avBlocks, d_numEigs * d_numBlocks, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // avfac already divides the number of samples per node by the total number -- so this is actually the average!
}

Discretizer::~Discretizer()
{
  delete[] d_hist;
  delete[] d_cum;
  delete[] d_block;
  delete[] d_avBlocks;
}
