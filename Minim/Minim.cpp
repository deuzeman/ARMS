#include <Minim.h>

double Minim::chiSq(Eigen::ArrayXd const &pars, size_t const iter)
{
  d_generator.calculate(pars, iter);
  Eigen::ArrayXXd dataRats = d_data.ratios(1000);
  double result = ((dataRats.row(0) - d_generator.ratios()) / dataRats.row(1)).square().sum(); // NOTE d_generator will always be symmetric, data won't be necessarily. Introduce that flexibility ?!
  return result;
}


double Minim::brent(Eigen::Array< double, 1, Eigen::Dynamic > const &center, Eigen::Array< double, 1, Eigen::Dynamic > const &dir, Eigen::ArrayXXd const &bounds, size_t const rmIters, double const tol)
{
  double const cgold = 0.3819660;

  Eigen::ArrayXXd limits(bounds);
  
  limits.row(0) = bounds.row(0) - center;
  limits.row(1) = bounds.row(1) - center; // Can be better, perhaps -- replicate somehow.

  double a = 0.0;
  double b = 0.0;
  double x = 0.0;

  // Bounds need to be calculated where the direction is non-zero.
  for (size_t col = 0; col < limits.cols(); ++col)
    for (size_t row = 0; row < limits.rows(); ++row)
    {
      if (dir(col) == 0.0)
	continue;
      double tmp = limits(row, col) / dir(col);
      if (tmp < 0)
        a = (tmp < a) ? tmp : a;
      else
	b = (tmp > b) ? tmp : b;
    }

  double u = 0.0;
  double w = 0.0;
  double v = 0.0;

  double fx = chiSq(center, rmIters);
  double fu = fx;
  double fv = fx;
  double fw = fx;

  double d = 0.0;
  double e = 0.0;

  for (size_t iter = 0; iter < 100; ++iter)
  {
    std::cerr << "[DEBUG -- BRENT] a, x, b, fx, fu: " << a << ' ' << x << ' ' << b << ' ' << fx << ' ' << fu << std::endl;

    double xm = 0.5 * (a + b);
    e = (x > xm) ? (a - x) : (b - x);
    if (std::abs(e) <= (tol))
      return x;

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
      u = x + d;
      if ((u - a) < (2 * tol) || (b - u) < (2 * tol))
	d = (xm <= x) ? -tol : tol;
    }

    u = (std::abs(d) >= tol) ? (x + d) : x + (d > 0 ? tol : -tol);
    fu = chiSq(center + u * dir, rmIters);

    if (fu <= fx)
    {
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
    for (size_t idx = 0; idx < n - 1; ++idx)
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