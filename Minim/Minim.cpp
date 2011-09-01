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

  std::cerr << "[DEBUG -- CHISQ ] " << padPars.transpose() << std::endl;
  
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
  size_t const maxIters = 400000;
  size_t const nBoot = 5000;
  
  Eigen::ArrayXXd limits(bounds);

  limits.row(0) = bounds.row(0) - center;
  limits.row(1) = bounds.row(1) - center; // Can be better, perhaps -- replicate somehow.

  std::cerr << "[DEBUG -- BRENT] Called Brent routine on:" << std::endl;
  std::cerr << "[DEBUG -- BRENT] " << center << std::endl;
  std::cerr << "[DEBUG -- BRENT] " << dir << std::endl;

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

  double fx = chiSq(center, rmIters, &ex, nBoot);
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
    fu = chiSq(center + u * dir, rmIters, &eu, nBoot);
    std::cerr << "[DEBUG -- BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;

    // Sufficient precision?
    std::cerr << "[DEBUG -- BRENT] Difference: " << fu - fx << " equiv " << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << std::endl;
    while (std::abs(fu - fx) < (3 * std::sqrt(eu * eu + ex * ex)))
    {
      std::cerr << "[DEBUG -- BRENT] Insufficient separation, additional precision needed." << std::endl;
      double fac = 1.33 * (3 * std::sqrt(eu * eu + ex * ex)) / std::abs(fu - fx);
      size_t newIters = static_cast< size_t >(fac * fac * static_cast< double >(rmIters));
      std::cerr << "[DEBUG -- BRENT] New value for rmIters: " << newIters << std::endl;
      if (newIters > maxIters)
      {
        std::cerr << "[DEBUG -- BRENT] This value is too high, changing to finding phi bounds." << std::endl;
        std::cerr << "[DEBUG -- BRENT] a = " << a << ", b = " << b << ", x = " << x << ", u = " << u << "." << std::endl;

        switchGen();
        fx = chiSq(center + x * dir, maxIters - rmIters, &ex, nBoot, true);
        std::cerr << "[DEBUG -- BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;

        switchGen();
        fu = chiSq(center + u * dir, maxIters - rmIters, &eu, nBoot, true);
        std::cerr << "[DEBUG -- BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
        std::cerr << "[DEBUG -- BRENT] Difference: " << fu - fx << " equiv " << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << std::endl;
        if (std::abs(fu - fx) > (3 * std::sqrt(eu * eu + ex * ex))) // Who knows? We may be out of the woods already...
          break;

        // First check that we're not accidentally on two sides of a proper minimum!
        double fb = 0.0;
        double eb = 0.0;
        double xb = 0.0;

        if (std::abs(x - u) > 1.5 * tol)
        {
          std::cerr << "[DEBUG -- BRENT] Checking for an internal minimum..." << std::endl;
          xb = phi(x, u, maxIters, center, dir, fx, ex, tol, &fb, &eb);
          if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex)))) // Best case scenario!
          {
            std::cerr << "[DEBUG -- BRENT] Phi bracketing produced a new minimum at " << xb << std::endl;
            // We'll pretend this was the current value of u and continue as if nothing happened.
            fu = fb;
            eu = eb;
            u  = xb;
            break;
          }
        }
        else
          std::cerr << "[DEBUG -- BRENT] No room for an internal minimum." << std::endl;

        std::cerr << "[DEBUG -- BRENT] Checking from x to the bound." << std::endl;
        double edge = (u > x) ? a : b; // Start in the opposite direction, because we haven't scanned that yet...

        xb = phi(x, edge, maxIters, center, dir, fx, ex, tol, &fb, &eb);

        if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex)))) // Best case scenario!
        {
          std::cerr << "[DEBUG -- BRENT] Phi bracketing produced a new minimum at " << xb << std::endl;
          // We'll pretend this was the current value of u and continue as if nothing happened.
          fu = fb;
          eu = eb;
          u  = xb;
          break;
        }

        // If we're here, we've either reached the edge (and this is now xb) or we've found some non-trivial bound.
        // Either way, we know x isn't well determined up to xb...

        std::cerr << "[DEBUG -- BRENT] Checking from u to the other bound." << std::endl;
        edge = (u > x) ? b : a; // Do the other side now, starting at u since we already know it's no different from x.
        double ub = phi(u, edge, maxIters, center, dir, fx, ex, tol, &fb, &eb);

        if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex)))) // Best case scenario again!
        {
          std::cerr << "[DEBUG -- BRENT] Phi bracketing produced a new minimum at " << xb << std::endl;
          // We'll pretend this was the current value of u and continue as if nothing happened.
          fu = fb;
          eu = eb;
          u  = ub;
          break;
        }

        // That pretty much ends our routine without a true minimum, so we'll return an appropriate value here.
        // Since it seems our starting point was as good as a minimum, let's just return it. This should minimize oscillations
        // in a multidimensional parameter space where all values are rougly equivalent.
        std::cerr << "[DEBUG -- BRENT] Found an acceptable region between " << xb << " and " << ub << std::endl;
        double ret = 0.0;

        if ((x < tol) && (x > -tol))
          std::cerr << "[DEBUG -- BRENT] This is compatible with the lack of effect seen so far, so return " << ret << " to stabilize." << std::endl;
        else
        {
          ret = round((xb + ub) / 2.0, tol);
          std::cerr << "[DEBUG -- BRENT] We have indications that the result is unequal zero anyway, so return the midpoint " << ret << '.' << std::endl;
        }

        return ret;

      }

      switchGen();
      fx = chiSq(center, newIters - rmIters, &ex, nBoot, true);
      std::cerr << "[DEBUG -- BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;

      switchGen();
      fu = chiSq(center, newIters - rmIters, &eu, nBoot, true);
      std::cerr << "[DEBUG -- BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
      std::cerr << "[DEBUG -- BRENT] Difference: " << fu - fx << " equiv " << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << std::endl;

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

void Minim::powell(Eigen::ArrayXd &start, Eigen::ArrayXXd &bounds, size_t const rmIters, size_t const powIters, double const tol)
{
  size_t const n = start.size();

  // Start with the parameters as principal directions
  Eigen::MatrixXd dirs = Eigen::MatrixXd::Identity(n, n);
  // Store succesive steps in P
  Eigen::MatrixXd pars = Eigen::MatrixXd::Zero(n, n + 1);

  round(start, tol);
  round(bounds, tol);
  round(dirs, tol);
  
  pars.col(0) = start;

  for (size_t iter = 0; iter < powIters; ++iter)
  {
    for (size_t idx = 0; idx < n; ++idx)
    {
      std::cerr << "[DEBUG -- POWELL] Before minimization:   " << pars.col(idx).transpose() << std::endl;
      std::cerr << "[DEBUG -- POWELL] Calling on direction:  " << dirs.col(idx).transpose() << std::endl;
      pars.col(idx + 1) = pars.col(idx) + brent(pars.col(idx), dirs.col(idx), bounds, rmIters, tol) * dirs.col(idx);
      round(pars, tol);
      std::cerr << "[DEBUG -- POWELL] After minimization:    " << pars.col(idx + 1).transpose()  << std::endl;
    }

    for (size_t idx = 0; idx < n - 2; ++idx)
      dirs.col(idx).swap(dirs.col(idx + 1));

    dirs.col(n - 1) = pars.col(n) - pars.col(0);
    
    Eigen::MatrixXd newP0(pars.col(n));
    if (dirs.col(n - 1).squaredNorm() > tol)
    {  
      dirs.col(n - 1).normalize(); // We work on absolute precision, so this matters for the minimization!
      round(dirs, tol);
      
      std::cerr << "[DEBUG -- POWELL] Before minimization:         " << pars.col(n).transpose() << std::endl;
      std::cerr << "[DEBUG -- POWELL] Calling on extra direction:  " << dirs.col(n - 1).transpose() << std::endl;
      newP0 += brent(pars.col(n), dirs.col(n - 1), bounds, rmIters, tol) * dirs.col(n - 1);
      round(newP0, tol);
      std::cerr << "[DEBUG -- POWELL] After minimization:          " << newP0.transpose() << std::endl;
    }

    if (((newP0 - pars.col(0)).squaredNorm()) < tol)
    {
      d_res = newP0;
      std::cerr << "[DEBUG -- POWELL] Minimization complete!" << std::endl; 
      return;
    }

    pars.col(0) = newP0; // Otherwise we just start from 0 again... :(

    // Reorthogonalize to avoid linear dependency build-up, pace Brent & NR.
    std::cerr << "[DEBUG -- POWELL] Before reorthogonalization:\n" << dirs << std::endl;
    Eigen::JacobiSVD< Eigen::MatrixXd > svd(dirs, Eigen::ComputeFullU);
    dirs = svd.matrixU();
    round(dirs, tol);
    std::cerr << "[DEBUG -- POWELL] After reorthogonalization:\n" << dirs << std::endl;
  }
  std::cerr << "[POWELL] Warning: Iteration count exceeded!" << std::endl;
}


double Minim::phi(double const start, double const edge, size_t const iter, Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir,
                  double const value, double const error, double const tol, double * const bval, double * const berr)
{
  double const cgold = 0.3819660;
  size_t const nBoot = 5000;

  std::cerr << "[DEBUG -- PHI] Entering bracketing routine to determine extent of plateau, between " << start << " and " << edge << '.' << std::endl;

  double attempt = round(start + 0.3819660 * (edge - start), tol);
  if (attempt < start + tol && attempt > start - tol) // i.e. equal to the starting point...
    attempt += (edge - start > 0) ? tol : -tol; // Now it's different at least

  Eigen::ArrayXd pars = center + attempt * dir;

  *bval = chiSq(pars, iter, berr, nBoot, false);

  if (std::abs(edge - attempt) < tol)
  {
    std::cerr << "[DEBUG -- PHI] Terminate on far edge condition." << std::endl;
    return attempt; // End of iterations anyway, because we've reached the edge...
  }

  if (std::abs(*bval - value) < (3 * std::sqrt(error * error + *berr * *berr))) // Still not there yet...
    return phi(attempt, edge, iter, center, dir, value, error, tol, bval, berr);

  if (*bval < value) // This is actually a new minimum!
  {
    std::cerr << "[DEBUG -- PHI] Terminate on new minimum." << std::endl;
    return attempt; // We're out of the slump, return to Brent.
  }

  // The only remaining (and most likely) case: we've found an increase in chi squared.

  if (std::abs(start - attempt) < tol) // No room for further improvement...
  {
    std::cerr << "[DEBUG -- PHI] Terminate on near edge condition." << std::endl;
    return attempt; // End of iterations anyway, because we've reached the edge...
  }

  // Room left, so go for it!
  return phi(start, attempt, iter, center, dir, value, error, tol, bval, berr);
}
