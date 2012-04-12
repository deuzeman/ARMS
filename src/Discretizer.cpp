#include <Discretizer.h>

#include <numeric>

Discretizer::Discretizer(double const *breaks, size_t numBreaks, size_t numEigs, size_t blocks)
 : d_numEigs(numEigs), d_levels(numBreaks + 1), d_numBlocks(blocks), d_breaks(breaks), d_data(0), d_hist(0), d_cum(0), d_block(0)
{
   d_acc = new double[d_numBlocks];
   d_hist = new double[d_levels * d_numBlocks];
   d_cum = new double[d_numEigs * d_levels];
   d_block = new double[d_numEigs * d_levels * d_numBlocks];
   
   clear();
}

void Discretizer::clear()
{
  std::fill_n(d_cum, d_numEigs * d_levels, 0.0);
  std::fill_n(d_block, d_numEigs * d_levels * blocks, 0.0);
  
  d_sampTotal = 0;
}
 
void Discretizer::calculate(RanMat const &ranmat)
{
  d_numEigs = ranmat.numEigs();
  d_numBlocks = d_numBlocks;
  d_levels = numBreaks + 1; 
  
  double const *res = ranmat.result();
  d_data = new size_t[ranmat.numSamples() * d_numEigs];
  
  for (size_t eig = 0; eig < d_numEigs; ++eig)
  {
    size_t *pd = d_data + eig * ranmat.numSamples();
    double const *pb = breaks + eig * numBreaks;
    for (size_t samp = 0; samp < ranmat.numSamples(); ++samp)
      while ((pd[samp] < numBreaks) && res[samp] > pb[pd[samp]])
        ++pd[samp];
  }

  // Check if we are extending our data here.
  // If so, fix the normalization of the data here.
  if (d_sampTotal > 0)
  {
    double refac = static_cast< double >(d_sampTotal) / (d_sampTotal + ranmat.numSamples());
    // d_cum will be overwritten later on, so we don't need to bother rescaling
    std::fill_n(d_cum, d_levels * d_numEigs, 0.0);
    for (size_t idx = 0; idx < d_numEigs * d_levels * d_numBlocks; ++idx)
      d_block[idx] *= refac;
  }
  d_sampTotal += ranmat.numSamples();
  
  double const inc = 1.0 / d_sampTotal;
  
  size_t const bStep = d_numEigs * d_levels;

  for (size_t eig = 0; eig < d_numEigs; ++eig)
  {
    std::fill_n(d_hist, d_levels * d_numBlocks, 0.0);

    size_t *pd = d_data + eig * ranmat.numSamples();
    double *pc = d_cum + eig * d_levels;
    
    for (size_t samp = 0; samp < ranmat.numSamples(); ++samp)
      d_hist[(samp % d_numBlocks) * d_levels + pd[samp]] += inc;
    
    std::fill_n(d_acc, d_numBlocks, 0.0);
    for (size_t lev = 0; lev < d_levels; ++lev)
    {
      for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
      {
        d_acc[bIdx] += d_hist[bIdx * d_levels + lev];
        d_block[(bIdx * d_numEigs + eig) * d_levels + lev] += d_acc[bIdx];
        pc[lev] += d_block[(bIdx * d_numEigs + eig) * d_levels + lev];
      }
      std::cerr << "[DEBUG] Final value in cumulative distribution: " << pc[d_levels - 1] << " (should be 1)." << std::endl;
    }
  }
}
    
  
  
  
  
  for (size_t bIdx = 0; bIdx < d_numBlocks; ++bIdx)
  {
    for (size_t eig = 0; eig < d_numEigs; ++eig)
    {
      double *pbh = d_block + bIdx * bStep + eig * d_levels;
      double *ph = d_hist + eig * d_levels;
      
      std::partial_sum(pbh, pbh + eig * d_levels, ph);
      std::cerr << "[DEBUG] Final value in block cumulative distribution: " << ph[d_levels - 1] << std::endl;
    }
    std::copy(d_hist, d_hist + d_numEigs * d_levels, d_block + bIdx * bStep);
  }
  
  delete[] d_data;
  d_data = 0;
}

Discretizer::~Discretizer()
{
  delete[] d_hist;
  delete[] d_cum;
  delete[] d_block;
}
