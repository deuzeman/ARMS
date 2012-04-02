#include <Minim.h>

Minim::Minim(Data const &data, size_t const N, size_t const nu, size_t const kol)
  : d_kol(kol), d_data(data), d_hash_gen(new Eigen::ArrayXXi()), d_hash_sec(new Eigen::ArrayXXi())
{
  MPI_Comm_rank(MPI_COMM_WORLD, &d_mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_numnodes);

  int blocklens[2] = {2, 5};
  MPI_Aint indices[2] = {0, 2 * sizeof(int)};
  MPI_Datatype old_types[2] = {MPI_INT, MPI_DOUBLE};
  MPI_Type_struct(2, blocklens, indices, old_types, &MPI_CONTROL);
  MPI_Type_commit(&MPI_CONTROL);

  d_countBuffer = new unsigned long int[d_data.numCols() * d_data.numSamples()];
  d_countBackup = new unsigned long int[d_data.numCols() * d_data.numSamples()]; // Since we can't reduce in place!

  if (d_mpirank == 0)
    d_log.open("fit.log", std::ios::app);

  // Set the appropriate seeds
  srand(time(0));
  int seeds[3] = {rand(), rand(), rand()};
  MPI_Bcast(&seeds, 3, MPI_INT, 0, MPI_COMM_WORLD);

  seeds[0] += d_mpirank;
  d_generator = new RanMat(N, nu, data.minEv(), data.maxEv(), 0, seeds[0]);

  seeds[1] += d_mpirank;
  d_secondary = new RanMat(N, nu, data.minEv(), data.maxEv(), 0, seeds[1]);

  d_sampler = new CRandomSFMT(seeds[2]); // We *want* identical sequences here
}

Minim::~Minim()
{
  d_log.close();

  MPI_Type_free(&MPI_CONTROL);

  delete d_generator;
  delete d_secondary;
  delete d_sampler;

  delete[] d_countBuffer;
  delete[] d_countBackup;
}

double Minim::chiSq(Eigen::ArrayXd const &pars, size_t const iter, double * const error, size_t const nBoot, bool const extend)
{
  d_log << "[CHISQ] " << pars.transpose() << std::endl;
  size_t ceilIters = iter + static_cast< size_t >(d_numnodes) - 1; 

  // Prepare the minions to support us here
  d_control.order[0] = extend ? 3 : 2;
  d_control.order[1] = static_cast< int >(ceilIters) / d_numnodes;
  d_control.params[0] = pars[0]; d_control.params[1] = pars[1]; d_control.params[2] = pars[2]; d_control.params[3] = pars[3]; d_control.params[4] = pars[4];
  MPI_Bcast(&d_control, 1, MPI_CONTROL, 0, MPI_COMM_WORLD);

  // The minions have scurried away to calculate their part. Let's take care of ours now.
  d_generator->calculate(pars, d_control.order[1], extend);
  buildHash(extend);
  clearCountBuffer();
  for (size_t idx = 0; idx < d_control.order[1]; ++idx)
    addSampleToCount(idx);
  condense();

  // Convert the simulation to a chi squared value now
  double result = 0.0;

  // Pull the branch point outside for performance reasons 
  if (d_kol)
  {
    for (size_t col = 0; col < d_data.numCols() ; ++col)
    {
      double pfrac = 0.0;
      for (size_t row = 0; row < d_data.numSamples(); ++row)
      {
        double dfrac = static_cast< double >(row) / static_cast< double >(d_data.numSamples());
        pfrac += static_cast< double >(d_countBuffer[col * d_data.numSamples() + row]) / static_cast< double >(ceilIters);
        result = std::max(result, std::abs(dfrac - pfrac));
      }
    }
    // No normalization needed (probably)
  }
  else
  {
    for (size_t col = 0; col < d_data.numCols() ; ++col)
    {
      double pfrac = 0.0;
      for (size_t row = 0; row < d_data.numSamples(); ++row)
      {
        double dfrac = static_cast< double >(row) / static_cast< double >(d_data.numSamples());
        pfrac += static_cast< double >(d_countBuffer[col * d_data.numSamples() + row]) / static_cast< double >(ceilIters);
        result += (dfrac - pfrac) * (dfrac - pfrac);
      }
    }
    result *= d_data.normalization();
  }

  // If an error pointer is provided, we need to bootstrap the result.
  if (error)
  {
    // Let the minions know they need to boostrap the RMT.
    // Sampling is done through taking the same random numbers from a synchronized state random number generator
    d_control.order[0] = 4; // BOOTSTRAP
    d_control.order[1] = nBoot;
    MPI_Bcast(&d_control, 1, MPI_CONTROL, 0, MPI_COMM_WORLD);

    size_t numSamples = d_hash_gen->rows() * d_numnodes;
    Eigen::ArrayXd bootHist(nBoot);

    for (size_t boot = 0; boot < nBoot; ++boot)
    {
      clearCountBuffer();
      for (size_t ctr = 0; ctr < numSamples; ++ctr)
      {
        size_t ranSample = static_cast< size_t >(d_sampler->IRandomX(0, numSamples - 1));
        if ((ranSample / d_hash_gen->rows()) == d_mpirank)
          addSampleToCount(ranSample % d_hash_gen->rows());
      }
      condense();

      bootHist[boot] = 0.0;
      // Pull the branch point outside for performance reasons 
      if (d_kol)
      {
        for (size_t col = 0; col < d_data.numCols() ; ++col)
        {
          double pfrac = 0.0;
          for (size_t row = 0; row < d_data.numSamples(); ++row)
          {
            double dfrac = static_cast< double >(row) / static_cast< double >(d_data.numSamples());
            pfrac += static_cast< double >(d_countBuffer[col * d_data.numSamples() + row]) / static_cast< double >(ceilIters);
            bootHist[boot] = std::max(bootHist[boot], std::abs(dfrac - pfrac));
          }
        }
        // No normalization needed (probably)
      }
      else
      {
        for (size_t col = 0; col < d_data.numCols() ; ++col)
        {
          double pfrac = 0.0;
          for (size_t row = 0; row < d_data.numSamples(); ++row)
          {
            double dfrac = static_cast< double >(row) / static_cast< double >(d_data.numSamples());
            pfrac += static_cast< double >(d_countBuffer[col * d_data.numSamples() + row]) / static_cast< double >(ceilIters);
            bootHist[boot] += (dfrac - pfrac) * (dfrac - pfrac);
          }
        }
        bootHist[boot] *= d_data.normalization();
      }
    }
    *error = std::sqrt(std::abs(bootHist.square().mean() - bootHist.mean() * bootHist.mean()));
  }
  return result;
}

void Minim::buildHash(bool const extend)
{
  Eigen::Array< double, Eigen::Dynamic, Eigen::Dynamic > const &samples = d_generator->samples();
  size_t offset = 0;

  if (extend)
  {
    offset = d_hash_gen->rows();
    d_hash_gen->conservativeResize(samples.rows(), samples.cols());
  }
  else
  {
    d_hash_gen->resize(samples.rows(), samples.cols());
  }

  Eigen::ArrayXXd const &data = d_data.cumulant();

  for (size_t col = 0; col < samples.cols(); ++col)
  {
    for (size_t row = offset; row < samples.rows(); ++row)
    {
      int stepper = 0; 
      while ((stepper < data.rows()) && (samples(row, col) > data(stepper, col)))
        ++stepper;
      (*d_hash_gen)(row, col) = stepper;
    }
  }
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

  double a = -1.0;
  double b = 1.0;
  double x = 0.0;
  double ex = 0.0;
  double eu = 0.0;

  double xcache = -1.1; // Will never occur...
  double fcache = 0.0;
  double ecache = 0.0;

  // Bounds need to be calculated where the direction is non-zero.
  for (size_t col = 0; col < limits.cols(); ++col)
  {
      if (dir(col) == 0.0)
        continue;

      size_t lowest = limits(0, col) < limits(1, col) ? 0 : 1;
      size_t highest = lowest ? 0 : 1;
      if (dir(col) > 0.0)
      {
        // Lower bound
        double tmp = limits(lowest, col) / dir(col);
        a = round((tmp > a) ? tmp : a, tol);

        // Higher bound
        tmp = limits(highest, col) / dir(col);
        b = round((tmp < b) ? tmp : b, tol);
      }
      else
      {
        // Lower bound
        double tmp = limits(highest, col) / dir(col);
        a = round((tmp > a) ? tmp : a, tol);

        // Higher bound
        tmp = limits(lowest, col) / dir(col);
        b = round((tmp < b) ? tmp : b, tol);
      }
  }

  // If the limits are too close, just return here.
  if (std::abs(b - a) < tol)
    return 0.0;

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
    if (u == xcache) // Does this ever happen? Don't think so, actually
    {
      d_log << "[BRENT] Cached results for u = " << u << " (x = " << x << ')' << std::endl;
      fu = fcache;
      eu = ecache;
    }
    else
    {
      d_log << "[BRENT] Calculating at u = " << u << " (x = " << x << ')' << std::endl;
      fu = chiSq(center + u * dir, rmIters, &eu, nBoot);
    }
    d_log << "[BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;

    // Sufficient precision?
    d_log << "[BRENT] Difference: " << fu - fx << " (" << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << " standard deviations)." << std::endl;
    bool phiMin = false;
    while (std::abs(fu - fx) < (3 * std::sqrt(eu * eu + ex * ex))) // Add condition to avoid rubbish points equivalent to zero...
    {
      d_log << "[BRENT] Insufficient separation, additional precision needed." << std::endl;
      double fac = 1.33 * (3 * std::sqrt(eu * eu + ex * ex)) / std::max(std::abs(fu - fx), 1E-10);
      size_t newIters = static_cast< size_t >(fac * fac * static_cast< double >(rmIters));
      if (newIters == 0)
        newIters = maxIters + 1;

      d_log << "[BRENT] New value for rmIters: " << newIters << std::endl;
      if (newIters > maxIters)
      {
        d_log << "[BRENT] This value is above the maximum, so we need to prepare for finding phi bounds." << std::endl;
        d_log << "[BRENT] a = " << a << ", b = " << b << ", x = " << x << ", u = " << u << "." << std::endl;

        switchGen();
        if (x == xcache)
        {
          d_log << "[BRENT] Cached results for x = " << x << ": " << (center + x * dir) << std::endl;
          fx = fcache;
          ex = ecache;
        }
        else
        {
          d_log << "[BRENT] Extending precision for x = " << x << ": " << (center + x * dir) << std::endl;
          fx = chiSq(center + x * dir, maxIters - rmIters, &ex, nBoot, true);
          xcache = x;
          fcache = fx;
          ecache = ex;
        }
        d_log << "[BRENT] Calculated fx +/- std as: " << fx << " +/- " << ex << std::endl;

        switchGen();
        if (u == xcache)
        {
          d_log << "[BRENT] Cached results for u = " << u << ": " << (center + u * dir) << std::endl;
          fu = fcache;
          eu = ecache;
        }
        else
        {
          d_log << "[BRENT] Extending precision for u = " << u << ": " << (center + u * dir) << std::endl;
          fu = chiSq(center + u * dir, maxIters - rmIters, &eu, nBoot, true);
        }

        d_log << "[BRENT] Calculated fu +/- std as: " << fu << " +/- " << eu << std::endl;
        d_log << "[BRENT] Difference: " << fu - fx << " (" << (fu - fx) / (std::sqrt(eu * eu + ex * ex)) << " standard deviations)." << std::endl;
        if (std::abs(fu - fx) > (3 * std::sqrt(eu * eu + ex * ex)) || (std::sqrt(eu * eu + ex * ex)) > (fu / 4.0)) // Who knows? We may be out of the woods already...
        {
          d_log << "[BRENT] Separation sufficient at highest precision, so cancel phi bound searching for now." << std::endl;
          if (fu < fx)
          {
            xcache = u;
            fcache = fu;
            ecache = eu;
          }
          break;
        }

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
      if (x == xcache)
      {
        d_log << "[BRENT] Cached results for x = " << x << ": " << (center + x * dir) << std::endl;
        fx = fcache;
        ex = ecache;
      }
      else
      {
        d_log << "[BRENT] Extending precision for x = " << x << ": " << (center + x * dir) << std::endl;
        fx = chiSq(center + x * dir, newIters - rmIters, &ex, nBoot, true);
      }
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

      // We'll have cached results for this new minimum, too
      xcache = u;
      fcache = fu;
      ecache = eu;
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
  d_log << "[POWELL] Start of calculation.\n"
        << "[POWELL] Starting from: " << start.transpose() << ".\n"
        << "[POWELL] Eigenvalue numbers: " << d_data.minEv() << " to " << d_data.maxEv() << ".\n"
        << "[POWELL] Number of samples: " << d_data.numSamples() << ".\n"
      << "[POWELL] Normalization is " << 1.0 / ((d_data.numCols()) * (0.5 - (Eigen::ArrayXd::LinSpaced(d_data.numSamples(), 1.0, static_cast< double >(d_data.numSamples())) / 
                                                                   static_cast< double >(d_data.numSamples()))).square().sum()) << '.' << std::endl;

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

  if (d_mpirank == 0)
  {
    d_control.order[0] = -1; // END_PROGRAM
    MPI_Bcast(&d_control, 1, MPI_CONTROL, 0, MPI_COMM_WORLD);
  }

  if (iter == (powIters - 1))
    d_log << "[POWELL] Warning: Iteration count exceeded!" << std::endl;

  d_log << "[POWELL] Final result: " << d_res.transpose() << std::endl;

  d_generator->calculate(d_res, maxIters);
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

inline void Minim::switchGen()
{
  if (d_mpirank == 0)
  {
    d_control.order[0] = 1; // SWITCH_GEN
    MPI_Bcast(&d_control, 1, MPI_CONTROL, 0, MPI_COMM_WORLD);
  }
  RanMat *tmp = d_generator;
  d_generator = d_secondary;
  d_secondary = tmp;

  Eigen::ArrayXXi *tmp_hash = d_hash_gen;
  d_hash_gen = d_hash_sec;
  d_hash_sec = tmp_hash;
}

void Minim::listen()
{
  if (d_mpirank == 0)
    return; // Should never be called on root anyway

  size_t nBoot;
  size_t iter;
  size_t numSamples;
  Eigen::ArrayXd pars(5); // To parse the parameters into

  while (1)
  {
    MPI_Bcast(&d_control, 1, MPI_CONTROL, 0, MPI_COMM_WORLD); // This is a passive call, await root's orders
    switch(d_control.order[0])
    {
      case -1: // END_CYCLE
        return; // End of program, break cycle
      case 1: // SWITCH_GEN
        switchGen();
        break;
      case 2: // CALC_NEW
        iter = static_cast< size_t >(d_control.order[1]);
        pars << d_control.params[0], d_control.params[1], d_control.params[2], d_control.params[3], d_control.params[4];
        d_generator->calculate(pars, d_control.order[1], false);
        buildHash(false);
        clearCountBuffer();
        for (size_t idx = 0; idx < d_control.order[1]; ++idx)
          addSampleToCount(idx);
        condense();
        break;
      case 3: // CALC_EXTEND
        iter = static_cast< size_t >(d_control.order[1]);
        pars << d_control.params[0], d_control.params[1], d_control.params[2], d_control.params[3], d_control.params[4];
        d_generator->calculate(pars, d_control.order[1], true);
        buildHash(true);
        clearCountBuffer();
        for (size_t idx = 0; idx < d_control.order[1]; ++idx)
          addSampleToCount(idx);
        condense();
        break;
      case 4: // BOOTSTRAP
        nBoot = d_control.order[1];
        numSamples = d_hash_gen->rows() * d_numnodes;
        for (size_t boot = 0; boot < nBoot; ++boot)
        {
          clearCountBuffer();
          for (size_t ctr = 0; ctr < numSamples; ++ctr)
          {
            size_t ranSample = static_cast< size_t >(d_sampler->IRandomX(0, numSamples - 1));
            if ((ranSample / d_hash_gen->rows()) == d_mpirank)
              addSampleToCount(ranSample % d_hash_gen->rows());
          }
          condense();
        }
        break;
      default: // Weird unknown value, should be a programming error. Better crash here.
        return;
    }
  }
}

void Minim::condense()
{
  MPI_Reduce(d_countBuffer, d_countBackup, static_cast< int >(d_data.numSamples() * d_data.numCols()), MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  unsigned long int *tmp = d_countBuffer;
  d_countBuffer = d_countBackup;
  d_countBackup = tmp;
}
