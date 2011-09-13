#include <Minim.h>

double Minim::chiSq(Eigen::ArrayXd const &pars, size_t const iter, double * const error, size_t const nBoot, bool const extend)
{
  if (d_cumulant)
    return chiSqCum(pars, iter, error, nBoot, extend);
  else
    return chiSqRat(pars, iter, error, nBoot, extend);
}

double Minim::chiSqRat(Eigen::ArrayXd const &pars, size_t const iter, double * const error, size_t const nBoot, bool const extend)
{
  // Potentially pad to add the flexibility for a 2D miminization.
  d_log << "[CHISQRAT] " << pars.transpose() << std::endl;

  d_generator->calculate(pars, iter, extend);
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

double Minim::chiSqCum(Eigen::ArrayXd const &pars, size_t const iter, double * const error, size_t const nBoot, bool const extend)
{
  d_log << "[CHISQCUM]" << pars.transpose() << std::endl;

  d_generator->calculate(pars, iter, extend);
  Eigen::ArrayXXd const &data = d_data.cumulant();
  Eigen::ArrayXXd const &samples = d_generator->samples();
  Eigen::ArrayXXd pred(samples.rows(), samples.cols());
  Eigen::Array< double, 1, Eigen::Dynamic > tempVec(samples.rows());

  for (size_t idx = 0; idx < pred.cols(); ++idx)
  {
    tempVec = samples.col(idx);
    std::sort(&tempVec[0], &tempVec[pred.rows() - 1] + 1);
    pred.col(idx) = tempVec;
  }

  double const normalization = 1.0 / data.cols() * (0.5 - (Eigen::ArrayXd::LinSpaced(data.cols(), 1.0, static_cast< double >(data.rows())) / static_cast< double >(data.rows() * data.rows()))).square().sum();

  double result = 0.0;

  for (size_t col = 0; col < data.cols(); ++col)
  {
    size_t stepper = 0;
    for (size_t row = 0; row < data.rows(); ++row)
    {
      while ((data(row, col) > pred(stepper, col)) && (stepper < (pred.rows() - 1)))
        ++stepper;

      double dfrac = static_cast< double >(row) / static_cast< double >(data.rows());
      double pfrac = static_cast< double >(stepper) / static_cast< double >(pred.rows());
      result += (dfrac - pfrac) * (dfrac - pfrac);
    }
    result *= normalization;
  }

  if (error != 0)
  {
    Eigen::ArrayXXd bootSamp(samples.rows(), samples.cols());
    Eigen::ArrayXd bootHist(nBoot);
    CRandomSFMT sampler(time(0));

    for (size_t bootCtr = 0; bootCtr < nBoot; ++bootCtr)
    {
      for (size_t sampCtr = 0; sampCtr < samples.rows(); ++sampCtr)
        bootSamp.row(sampCtr) = samples.row(sampler.IRandomX(0, samples.rows() - 1));

      for (size_t idx = 0; idx < pred.cols(); ++idx)
      {
        tempVec = bootSamp.col(idx);
        std::sort(&tempVec[0], &tempVec[bootSamp.rows() - 1] + 1);
        bootSamp.col(idx) = tempVec;
      }

      bootHist[bootCtr] = 0.0;
      for (size_t col = 0; col < data.cols(); ++col)
      {
        size_t stepper = 0;
        for (size_t row = 0; row < data.rows(); ++row)
        {
          while ((data(row, col) > bootSamp(stepper, col)) && (stepper < (bootSamp.rows() - 1)))
            ++stepper;

          double dfrac = static_cast< double >(row) / static_cast< double >(data.rows());
          double pfrac = static_cast< double >(stepper) / static_cast< double >(bootSamp.rows());

          bootHist[bootCtr] += (dfrac - pfrac) * (dfrac - pfrac);
        }
        bootHist[bootCtr] *= normalization;
      }
    }

    *error = std::sqrt(std::abs(bootHist.square().mean() - bootHist.mean() * bootHist.mean()));
  }

  return result;
}


double Minim::brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, size_t rmIters, size_t const maxIters, double const tol)
{
  double const cgold = 0.3819660;
  size_t const nBoot = 5000;

  Eigen::ArrayXXd limits(bounds);

  limits.row(0) = bounds.row(0) - center;
  limits.row(1) = bounds.row(1) - center; // Can be better, perhaps -- replicate somehow.

  d_log << "[BRENT] Called Brent routine on:" << std::endl;
  d_log << "[BRENT] " << center << std::endl;
  d_log << "[BRENT] " << dir << std::endl;

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
      double tmp = (limits(row, col) - center(col)) / dir(col);
      if (tmp < 0)
        a = round((tmp < a) ? tmp : a, tol);
      else
	b = round((tmp > b) ? tmp : b, tol);
    }

  // If the limits are too close, just return here.
  if (std::abs(b - a) < tol)
    return 0,0;

  d_log << "[BRENT] Limits determined as: [" << a << ", " << b << "]." << std::endl;

  double u = 0.0;
  double w = 0.0;
  double v = 0.0;

  double fx = chiSq(center, rmIters, &ex, nBoot);
  d_log << "[BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;
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
      d_log << "[BRENT] Minimization completed!" << std::endl;
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
    d_log << "[BRENT] Calculating at u = " << u << " (x = " << x << ')' << std::endl;
    fu = chiSq(center + u * dir, rmIters, &eu, nBoot);
    d_log << "[BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;

    // Sufficient precision?
    d_log << "[BRENT] Difference: " << fu - fx << " (" << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << " standard deviations)." << std::endl;
    bool phiMin = false;
    while (std::abs(fu - fx) < (3 * std::sqrt(eu * eu + ex * ex))) // Add condition to avoid rubbish points equivalent to zero...
    {
      d_log << "[BRENT] Insufficient separation, additional precision needed." << std::endl;
      double fac = 1.33 * (3 * std::sqrt(eu * eu + ex * ex)) / std::max(std::abs(fu - fx), tol);
      size_t newIters = static_cast< size_t >(fac * fac * static_cast< double >(rmIters));
      if (newIters == 0)
        newIters = maxIters + 1;

      d_log << "[BRENT] New value for rmIters: " << newIters << std::endl;
      if (newIters > maxIters)
      {
        d_log << "[BRENT] This value is too high, changing to finding phi bounds." << std::endl;
        d_log << "[BRENT] a = " << a << ", b = " << b << ", x = " << x << ", u = " << u << "." << std::endl;

        switchGen();
        d_log << "[BRENT] Extending precision for x = " << x << ": " << (center + x * dir) << std::endl;
        fx = chiSq(center + x * dir, maxIters - rmIters, &ex, nBoot, true);
        d_log << "[BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;

        switchGen();
        d_log << "[BRENT] Extending precision for u = " << u << ": " << (center + u * dir) << std::endl;
        fu = chiSq(center + u * dir, maxIters - rmIters, &eu, nBoot, true);
        d_log << "[BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
        d_log << "[BRENT] Difference: " << fu - fx << " (" << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << " standard deviations)." << std::endl;
        if (std::abs(fu - fx) > (3 * std::sqrt(eu * eu + ex * ex)) || (std::sqrt(eu * eu + ex * ex)) > (fu / 4.0)) // Who knows? We may be out of the woods already...
          break;

        // First check that we're not accidentally on two sides of a proper minimum!
        double fb = 0.0;
        double eb = 0.0;
        double xb = 0.0;

        if (std::abs(x - u) > 1.5 * tol)
        {
          d_log << "[BRENT] Checking for an internal minimum..." << std::endl;
          xb = phi(x, u, maxIters, center, dir, fx, ex, tol, &fb, &eb);
          if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex))) && (std::sqrt(eb * eb + ex * ex)) < (fb / 4.0)) // Best case scenario!
          {
            d_log << "[BRENT] Phi bracketing produced a new minimum at " << xb << std::endl;

            // We'll pretend this was the current value of u and continue as if nothing happened.
	    // Except that we know that the minimum will be between u and x!
            a  = std::min(x, u);
            b  = std::max(x, u);
            fu = fb;
            eu = eb;
            u  = xb;
	    phiMin = true;
            break;
          }
        }
        else
          d_log << "[BRENT] No room for an internal minimum." << std::endl;

        d_log << "[BRENT] Checking from x to the nearest bound." << std::endl;
        double edge = (u > x) ? a : b; // Start in the opposite direction, because we haven't scanned that yet...

        xb = phi(x, edge, maxIters, center, dir, fx, ex, tol, &fb, &eb);

        if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex))) && (std::sqrt(eb * eb + ex * ex)) < (fb / 4.0)) // Best case scenario!
        {
          d_log << "[BRENT] Phi bracketing produced a new minimum at " << xb << std::endl;
          // We'll pretend this was the current value of u and continue as if nothing happened.
          // Except that we know that the minimum will be between the edge and x!
          a  = std::min(x, edge);
          b  = std::max(x, edge);
          fu = fb;
          eu = eb;
          u  = xb;
	  phiMin = true;
          break;
        }

        // If we're here, we've either reached the edge (and this is now xb) or we've found some non-trivial bound.
        // Either way, we know x isn't well determined up to xb...

        d_log << "[BRENT] Checking from u to the other bound." << std::endl;
        edge = (u > x) ? b : a; // Do the other side now, starting at u since we already know it's no different from x.
        double ub = phi(u, edge, maxIters, center, dir, fx, ex, tol, &fb, &eb);

        if ((fb < fx) && (std::abs(fb - fx) > (3 * std::sqrt(eb * eb + ex * ex))) && (std::sqrt(eb * eb + ex * ex)) < (fb / 4.0)) // Best case scenario again!
        {
          d_log << "[BRENT] Phi bracketing produced a new minimum at " << ub << std::endl;
          // We'll pretend this was the current value of u and continue as if nothing happened.
          // Except that we know that the minimum will be between the edge and u!
          a  = std::min(u, edge);
          b  = std::max(u, edge);
          fu = fb;
          eu = eb;
          u  = ub;
	  phiMin = true;
          break;
        }

        // That pretty much ends our routine without a true minimum, so we'll return an appropriate value here.
        // Since it seems our starting point was as good as a minimum, let's just return it. This should minimize oscillations
        // in a multidimensional parameter space where all values are rougly equivalent.
        d_log << "[BRENT] Found an acceptable region between " << xb << " and " << ub << std::endl;
        double ret = 0.0;

        if ((x < tol) && (x > -tol))
          d_log << "[BRENT] This is compatible with the lack of effect seen so far, so return " << ret << " to stabilize." << std::endl;
        else
        {
          ret = round((xb + ub) / 2.0, tol);
          d_log << "[BRENT] We have indications that the result is unequal zero anyway, so return the midpoint " << ret << '.' << std::endl;
        }

        return ret;
      }

      switchGen();
      d_log << "[BRENT] Extending precision for x = " << x << ": " << (center + x * dir) << std::endl;
      fx = chiSq(center + x * dir, newIters - rmIters, &ex, nBoot, true);
      d_log << "[BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;

      switchGen();
      d_log << "[BRENT] Extending precision for u = " << u << ": " << (center + u * dir) << std::endl;
      fu = chiSq(center + u * dir, newIters - rmIters, &eu, nBoot, true);
      d_log << "[BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
      d_log << "[BRENT] Difference: " << fu - fx << " (" << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << " standard deviations)." << std::endl;

      rmIters = newIters;
    }

    if (phiMin) // To avoid wasting time on larger bounds when the phi routine has already carved out an area
    {
      switchGen();
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    }

    if (fu <= fx)
    {
      d_log << "[BRENT] Found a new minimum!" << std::endl;
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
  d_log << "[BRENT] Warning: Iteration count exceeded!" << std::endl;
  return x;
}

void Minim::powell(Eigen::ArrayXd &start, Eigen::ArrayXXd &bounds, size_t const rmIters, size_t const maxIters, size_t const powIters, double const tol)
{
  size_t const n = start.size();

  // Start with the parameters as principal directions
  Eigen::MatrixXd dirs = Eigen::MatrixXd::Identity(n, n);
  // We want to fit the scale first for cumulative fits
  if (n % 2)
  {
    dirs.col(0).swap(dirs.col(1));
    dirs.col(0).swap(dirs.col(2));
    if (n == 5)
    {
      dirs.col(0).swap(dirs.col(3));
      dirs.col(0).swap(dirs.col(4));
    }
  }

  // Store succesive steps in P
  Eigen::MatrixXd pars = Eigen::MatrixXd::Zero(n, n + 1);

  round(start, tol);
  round(bounds, tol);
  round(dirs, tol);

  pars.col(0) = start;

  size_t iter = 0;

  for (; iter < powIters; ++iter)
  {
    for (size_t idx = 0; idx < n; ++idx)
    {
      d_log << "[POWELL] Before minimization:   " << pars.col(idx).transpose() << std::endl;
      d_log << "[POWELL] Calling on direction:  " << dirs.col(idx).transpose() << std::endl;
      pars.col(idx + 1) = pars.col(idx) + brent(pars.col(idx), dirs.col(idx), bounds, rmIters, maxIters, tol) * dirs.col(idx);
      round(pars, tol);
      d_log << "[POWELL] After minimization:    " << pars.col(idx + 1).transpose()  << std::endl;
    }

    for (size_t idx = 0; idx < n - 2; ++idx)
      dirs.col(idx).swap(dirs.col(idx + 1));

    dirs.col(n - 1) = pars.col(n) - pars.col(0);

    Eigen::MatrixXd newP0(pars.col(n));
    if (dirs.col(n - 1).squaredNorm() > (tol * tol))
    {  
      dirs.col(n - 1).normalize(); // We work on absolute precision, so this matters for the minimization!
      round(dirs, tol);

      d_log << "[POWELL] Before minimization:         " << pars.col(n).transpose() << std::endl;
      d_log << "[POWELL] Calling on extra direction:  " << dirs.col(n - 1).transpose() << std::endl;
      newP0 += brent(pars.col(n), dirs.col(n - 1), bounds, rmIters, maxIters, tol) * dirs.col(n - 1);
      round(newP0, tol);
      d_log << "[POWELL] After minimization:          " << newP0.transpose() << std::endl;
    }

    if (((newP0 - pars.col(0)).squaredNorm()) < (tol * tol))
    {
      d_res = newP0;
      d_log << "[POWELL] Minimization complete!" << std::endl; 
      break;
    }

    pars.col(0) = newP0; // Otherwise we just start from 0 again... :(

    // Reorthogonalize to avoid linear dependency build-up, pace Brent & NR.
    d_log << "[POWELL] Before reorthogonalization:\n" << dirs << std::endl;
    Eigen::JacobiSVD< Eigen::MatrixXd > svd(dirs, Eigen::ComputeFullU);
    dirs = svd.matrixU();
    round(dirs, tol);
    d_log << "[POWELL] After reorthogonalization:\n" << dirs << std::endl;
  }
  
  if (iter == (powIters - 1))
    d_log << "[POWELL] Warning: Iteration count exceeded!" << std::endl;
  
  d_log << "[POWELL] Final result: " << d_res.transpose() ;

  d_generator->calculate(d_res, maxIters);

  std::ofstream rstream("fit_result.dat", std::ofstream::trunc);

  rstream << "# Params: " << d_res.transpose() << std::endl;

  for (int x = d_nEig_min; x < d_nEig_max; ++x)
  {
    int colNum = (x < 0) ? x : x + 1;
    std::ostringstream colLab;
    if (x < 0)
      colLab << "\"EV.m" << -colNum << '\"';
    else
      colLab << "\"EV.p" << colNum << '\"';
   rstream << std::setw(15) << colLab.str();
  }
  rstream << std::endl;

  for (int k = 0; k < maxIters; ++k)
  {
    std::ostringstream lineLab;
    lineLab << '\"' << (k + 1) << '\"';
    rstream << std::setw(10) << lineLab.str();

    for (int x = 0; x < d_nEig_max - d_nEig_min; ++x)
      rstream << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(8) << d_generator->result(k, x);

    rstream << std::endl;
  }
  rstream.close();
}


double Minim::phi(double const start, double const edge, size_t const iter, Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir,
                  double const value, double const error, double const tol, double * const bval, double * const berr)
{
  double const cgold = 0.3819660;
  size_t const nBoot = 5000;

  d_log << "[PHI] Entering bracketing routine to determine extent of plateau, between " << start << " and " << edge << '.' << std::endl;

  double attempt = round(start + 0.3819660 * (edge - start), tol);
  if (std::abs(attempt - start) < tol) // i.e. equal to the starting point...
    attempt += (edge - start > 0) ? tol : -tol; // Now it's different at least

  Eigen::ArrayXd pars = center + attempt * dir;
  d_log << "[PHI] Calculating with maximum precision for b = " << attempt << ": " << pars.transpose() << '.' <<  std::endl;

  *bval = chiSq(pars, iter, berr, nBoot, false);
  d_log << "[PHI] Calculated fb +/- std as: " << *bval << " +/- " << *berr << std::endl;
  d_log << "[PHI] Difference of " << (*bval - value) << " (" << std::abs(*bval - value) / std::sqrt(error * error + *berr * *berr) << " std)." << std::endl;

  if (std::abs(edge - attempt) < tol)
  {
    d_log << "[PHI] Terminate on far edge condition." << std::endl;
    return attempt; // End of iterations anyway, because we've reached the edge...
  }

  if (std::abs(*bval - value) < (3 * std::sqrt(error * error + *berr * *berr)) && (std::sqrt(error * error + *berr * *berr) < 0.5 * *bval)) // Still not there yet...
    return phi(attempt, edge, iter, center, dir, value, error, tol, bval, berr);

  if (*bval < value && (std::sqrt(error * error + *berr * *berr) < 0.5 * *bval)) // This is actually a new minimum!
  {
    d_log << "[PHI] Terminate on new minimum." << std::endl;
    return attempt; // We're out of the slump, return to Brent.
  }

  // The only remaining (and most likely) case: we've found an increase in chi squared.
  // Finding a relatively dominant error counts as such as well.


  // We've just calculated at attempt and found an increase in chi-squared.
  // No need to include it again as the edge -- it's out. But I do observe it as being calculated again...
  // To avoid wasting time, decrease the edge by one.
  // NOTE I believe we can assume start has been calculated already
  attempt += (attempt - start > 0) ? -tol : tol; // Now it's different at least

  if (std::abs(start - attempt) < tol) // No room for further improvement...
  {
    d_log << "[PHI] Terminate on near edge condition." << std::endl;
    return attempt; // End of iterations anyway, because we've reached the edge...
  }

  // Room left, so go for it!
  return phi(start, attempt, iter, center, dir, value, error, tol, bval, berr);
}
