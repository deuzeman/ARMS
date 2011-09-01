#include <Minim.h>

double Minim::chiSq(Eigen::ArrayXd const &pars, size_t const iter, double * const error, size_t const nBoot, bool const extend)
{
  // Potentially pad to add the flexibility for a 2D miminization
  Eigen::ArrayXd padPars = Eigen::ArrayXd::Zero(4);
  if (pars.size() == 2)
  {
    padPars[0] = pars[0];
    padPars[3] = pars[1];
  }
  else
  {
    padPars = pars;
  }
  
    
  d_generator->calculate(padPars, iter, extend);
  Eigen::ArrayXXd const &dataRats = d_data.ratios(1000);
  double result = ((dataRats.row(0) - d_generator->ratios().row(0)) / dataRats.row(1)).square().sum();

  if (error != 0) // Request for error calculation implied
  {
    Eigen::ArrayXXd const &samples = d_generator->samples();
    Eigen::ArrayXXd ratFull(samples.rows(), samples.cols() * (samples.cols() - 1) / 2);

    size_t ctr = 0;
    for (size_t num = 0; num < samples.cols() - 1; ++num)
      for (size_t den = num + 1; den < samples.cols(); ++den)
        ratFull.col(ctr++) = (samples.col(num) / samples.col(den));

    Eigen::ArrayXXd bootSamp(ratFull.rows(), ratFull.cols());
    Eigen::ArrayXd bootHist(nBoot);
    CRandomSFMT sampler(time(0));

    for (size_t bootCtr = 0; bootCtr < nBoot; ++bootCtr)
    {
      for (size_t sampCtr = 0; sampCtr < ratFull.rows(); ++sampCtr)
        bootSamp.row(sampCtr) = ratFull.row(sampler.IRandomX(0, ratFull.rows() - 1));
      bootHist[bootCtr] = ((dataRats.row(0) - (bootSamp.colwise().mean())) / dataRats.row(1)).square().sum();
    }
    *error = std::sqrt(bootHist.square().mean() - bootHist.mean() * bootHist.mean());
  }

  return result;
}


double Minim::brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, size_t rmIters, double const tol)
{
  double const cgold = 0.3819660;

  Eigen::ArrayXXd limits(bounds);

  limits.row(0) = bounds.row(0) - center;
  limits.row(1) = bounds.row(1) - center; // Can be better, perhaps -- replicate somehow.

  double a = 0.0;
  double b = 0.0;
  double x = 0.0;
  double ex = 0.0;
  double eu = 0.0;

  // Bounds need to be calculated where the direction is non-zero.
  for (size_t col = 0; col < limits.cols(); ++col)
    for (size_t row = 0; row < limits.rows(); ++row)
    {
      if (dir(col) == 0.0)
	continue;
      double tmp = limits(row, col) / dir(col);
      if (tmp < 0)
        a = round((tmp < a) ? tmp : a, tol);
      else
	b = round((tmp > b) ? tmp : b, tol);
    }

  double u = 0.0;
  double w = 0.0;
  double v = 0.0;

  double fx = chiSq(center, rmIters, &ex);
  std::cerr << "[DEBUG -- BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;
  switchGen();
  
  double fu = fx;
  double fv = fx;
  double fw = fx;

  double d = 0.0;
  double e = 0.0;

  for (size_t iter = 0; iter < 100; ++iter)
  {
    double xm = round(0.5 * (a + b), tol);
    e = (x > xm) ? (a - x) : (b - x);
    if (std::abs(e) < 1.5 * tol)
    {
      std::cerr << "[DEBUG -- BRENT] Minimization completed!" << std::endl;
      return x;
    }

    double r = (x - w) * (fx - fv);
    double q = (x - v) * (fx - fw);
    double p = (x - v) * q - (x - w) * r;
    q = 2.0 * (q - r);
    if (q > 0.0)
      p = -p;
    q = std::abs(q);
    double e2 = e;
    e = d;
    if (std::abs(p) >= std::abs(0.5 * q * e2) || p <= q * (a - x) || p >= q * (b - x))
    {
      e = (x > xm) ? (a - x) : (b - x);
      d = cgold * e;
    }
    else
    {
      d = p / q;
      u = round(x + d, tol);
      if ((u - a) < (2 * tol) || (b - u) < (2 * tol))
	d = (xm <= x) ? -tol : tol;
    }

    u = round((std::abs(d) >= tol) ? (x + d) : x + (d > 0 ? tol : -tol), tol);
    std::cerr << "[DEBUG -- BRENT] Calculating at u = " << u << " (x = " << x << ')' << std::endl;
    fu = chiSq(center + u * dir, rmIters, &eu);
    std::cerr << "[DEBUG -- BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;

    // Sufficient precision?
    std::cerr << "[DEBUG -- BRENT] Difference: " << fu - fx << " equiv " << (fu - fx) / (3 * std::sqrt(eu * eu + ex * ex)) << std::endl;
    while (std::abs(fu - fx) < (3 * std::sqrt(eu * eu + ex * ex)))
    {
      std::cerr << "[DEBUG -- BRENT] Not sufficient separation, additional precision needed." << std::endl;
      double fac = 1.33 * (3 * std::sqrt(eu * eu + ex * ex)) / std::abs(fu - fx);
      size_t newIters = static_cast< size_t >(fac * fac * static_cast< double >(rmIters));
      std::cerr << "[DEBUG -- BRENT] New value for rmIters: " << newIters << std::endl;
      if (newIters > 100000) // NOTE Hard-coded limit here...
      {
	std::cerr << "[DEBUG -- BRENT] Too many iterations! Reduce precision requirements..." << std::endl;
	exit(-1);
      }
      
      switchGen();
      fx = chiSq(center, newIters - rmIters, &ex, 1000, true);
      std::cerr << "[DEBUG -- BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;
      
      switchGen();
      fu = chiSq(center, newIters - rmIters, &eu, 1000, true);
      std::cerr << "[DEBUG -- BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
      std::cerr << "[DEBUG -- BRENT] Difference: " << fu - fx << " equiv " << (fu - fx) / (3 * std::sqrt(eu * eu + ex * ex)) << std::endl;

      rmIters = newIters;
    }
    
    if (fu <= fx)
    {
      std::cerr << "[DEBUG -- BRENT] Found a new minimum!" << std::endl;
      switchGen();
      if (u >= x)
        a = x;
      else
        b = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    }
    else
    {
      if (u < x)
        a = u;
      else
        b = u;
      if ((fu <= fw) || (w == x))
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else
        if ((fu <= fv) || (v == x) || (v == w))
        {
          v = u;
          fv = fu;
        }
    }
  }
  std::cerr << "[BRENT] Warning: Iteration count exceeded!" << std::endl;
  return x;
}

void Minim::powell(Eigen::VectorXd const &start, Eigen::ArrayXXd const &bounds, size_t const rmIters, size_t const powIters, double const tol)
{
  double const tiny = 1.0E-25;
  size_t const n = start.size();

  // Start with the parameters as principal directions
  Eigen::MatrixXd dirs = Eigen::MatrixXd::Identity(n, n);
  // Store succesive steps in P
  Eigen::MatrixXd pars = Eigen::MatrixXd::Zero(n, n + 1);

  pars.col(0) = start;

  for (size_t iter = 0; iter < powIters; ++iter)
  {
    for (size_t idx = 0; idx < n; ++idx)
    {
      std::cerr << "[DEBUG -- POWELL] Calling on direction:  " << dirs.col(idx).transpose() << std::endl;
      pars.col(idx + 1) = pars.col(idx) + brent(pars.col(idx), dirs.col(idx), bounds, rmIters, tol) * dirs.col(idx);
    }

    for (size_t idx = 0; idx < n - 2; ++idx)
      dirs.col(idx).swap(dirs.col(idx + 1));

    dirs.col(n - 1) = pars.col(n) - pars.col(0);
    Eigen::MatrixXd newP0 = pars.col(n) + brent(pars.col(n), dirs.col(n - 1), bounds, rmIters, tol) * dirs.col(n - 1);

    if (((newP0 - pars.col(0)).squaredNorm()) < tol)
    {
      d_res = newP0;
      return;
    }

    // Reorthogonalize to avoid linear dependency build-up, pace Brent & NR.
    Eigen::JacobiSVD< Eigen::MatrixXd > svd(dirs, Eigen::ComputeFullU);

    dirs = svd.matrixU();
  }
  std::cerr << "[POWELL] Warning: Iteration count exceeded!" << std::endl;
}