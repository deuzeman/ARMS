#include <iostream>
#include <cstdlib>
#include <cmath>
#include <numeric>

#include <mpi.h>

#include <Data.h>
#include <Params.h>
#include <RanMat.h>

double pks(double D, size_t N1, size_t N2)
{
  double Ne = (static_cast< double >(N1) * N2) / (N1 + N2);
  double arg = (std::sqrt(Ne) + 0.12 + 0.11 / std::sqrt(Ne)) * D;

  if (arg == 0)
    return 1.0;
  if (arg < 1.18)
  {
    double y = exp(-1.23370055013616983 / std::sqrt(arg));
    return 1.0 - (2.25675833419102515 * std::sqrt(-std::log(y)) * (y + std::pow(y, 9.0) + std::pow(y, 25.0) + std::pow(y, 49.0)));
  }
  else
  {
    double x = exp(-2 * std::sqrt(arg));
    return 2.0 * (x - std::pow(x, 4.0) + pow(x, 9.0));
  }
}

int main(int argc, char **argv)
{
  if (argc != 3)
  {
    std::cerr << "# RMGOF analysis code\n"
              << "Usage:\n" 
              << "  " << argv[0] << " [data_file] [args_input]\n"
	      << "Now quitting." << std::endl;
    return 1;
  }

  Data data(argv[1], 0, false);
  Params par(argv[2]);
  RanMat generator(par.N, par.nu, par.nEig_min, par.nEig_max, par.nDet, time(0), true);

  Eigen::ArrayXd rmtPars(5);
  rmtPars << par.m, par.a6, par.a7, par.a8, par.sigma;
  generator.calculate(rmtPars, par.iter, false);

  size_t ed = data.numCols() * data.numSamples();
  size_t ev = par.iter * (par.nEig_max - par.nEig_min);

  double *vd = data.flat();
  double *vr = generator.flat();
  
  size_t *breaks = new size_t[ed];

  size_t step = 0;
  for (size_t idx = 0; idx < ed; ++idx)
  {
    while (step < ev && (vr[step] < vd[idx])) ++step;
    breaks[idx] = step;
  }

  double *frac  = new double[ed];
  for (size_t idx = 0; idx < ed; ++idx)
    frac[idx] = static_cast< double >(breaks[idx]) / ev;

  double ks = 0.0;

  for (size_t idx = 0; idx < ed; ++idx)
    ks = std::max(ks, std::abs(frac[idx] - idx / static_cast< double >(ed)));

  std::cout << "# Using data from:                        " << argv[1] << std::endl;
  std::cout << "# Fit parameters:                         " << par.m << ' ' << par.a6 << ' ' << par.a7 << ' ' << par.a8 << ' ' << par.sigma << std::endl;
  std::cout << "# Kolmogorov-Smirnov deviation D:         " << ks << std::endl;
  std::cout << "# Probability of observing deviation > D: " << pks(ks, ev, ed) << std::endl;
  std::cout << "# " << std::endl;
 
  delete[] vd;
  delete[] vr;
  delete[] breaks;
  delete[] frac;

  ed = data.numSamples();
  ev = par.iter;
  double norm = 1 / static_cast< double >(ed);

  size_t **breaks_ar = new size_t*[data.numCols()];
  double **frac_ar = new double*[data.numCols()];
  double **vd_ar = new double*[data.numCols()];
  double **sd_low_ar = new double*[data.numCols()];
  double **sd_high_ar = new double*[data.numCols()];
  
  double *sd_low;
  double *sd_high;

  double *sd_work = new double[100 * ed];
  double *sd_work_vec = new double[ed];
  
  std::cout << "# ========================================" << std::endl;

  // Now analyze each eigenvalue separately

  for (size_t col = 0; col < data.numCols(); ++col)
  {
    vr = generator.flat(col);

    breaks_ar[col] = new size_t[ed];
    breaks = breaks_ar[col];

    frac_ar[col] = new double[ed];
    frac = frac_ar[col];

    vd_ar[col] = data.flat(col);
    vd = vd_ar[col];

    sd_low_ar[col] = new double[ed];
    sd_low = sd_low_ar[col];
    sd_high_ar[col] = new double[ed];
    sd_high = sd_high_ar[col];
    
    CRandomSFMT sampler(time(0));
    for (size_t idx = 0; idx < 100 * ed; ++idx)
      sd_work[idx] = vd[sampler.IRandomX(0, ed - 1)];
    
    for (size_t idx = 0; idx < 100; ++idx)
      std::sort(sd_work + idx * ed, sd_work + (idx + 1) * ed);
    
    for (size_t idx = 0; idx < ed; ++idx)
    {
      for (size_t row = 0; row < 100; ++row)
        sd_work_vec[row] = sd_work[idx + row * ed];
      std::sort(sd_work_vec, sd_work_vec + 100);
      sd_low[idx] = sd_work_vec[4];
      sd_high[idx] = sd_work_vec[95];
    }

    size_t step = 0;
    for (size_t idx = 0; idx < ed; ++idx)
    {
      while (step < ev && (vr[step] < vd[idx])) ++step;
      breaks[idx] = step;
    }

    for (size_t idx = 0; idx < ed; ++idx)
      frac[idx] = static_cast< double >(breaks[idx]) / ev;

    double ks = 0.0;

    for (size_t idx = 0; idx < ed; ++idx)
      ks = std::max(ks, std::abs(frac[idx] - norm * idx));

    std::cout << "# Peak " << col << ", D: " << ks << ", p-value: " << pks(ks, ev, ed) << std::endl;
  }
  
  std::cout << "# ========================================" << std::endl;
  std::cout << '#' << std::endl;
  std::cout << "   col  ev   sd_ev_min   sd_ev_max   cdf   fit" << std::endl;
  
  size_t counter = 0;
  for (size_t col = 0; col < data.numCols(); ++col)
  {
    vd = vd_ar[col];
    frac = frac_ar[col];
    sd_low = sd_low_ar[col];
    sd_high = sd_high_ar[col];
    for (size_t idx = 0; idx < ed; ++idx)
      std::cout << '\"' << ++counter << "\"  " << col << "  " << vd[idx] << "  " << sd_low[idx] << "  " << sd_high[idx] << "  " << norm * idx << "  " << frac[idx] << std::endl;
  }

  delete[] sd_work;
  delete[] sd_work_vec;

  for (size_t col = 0; col < data.numCols(); ++col)
  {
    delete[] vd_ar[col];
    delete[] breaks_ar[col];
    delete[] frac_ar[col];
    delete[] sd_low_ar[col];
    delete[] sd_high_ar[col];
  }
  delete[] vd_ar;
  delete[] breaks_ar;
  delete[] frac_ar;
  delete[] sd_low_ar;
  delete[] sd_high_ar;
  
  return 0;
}
