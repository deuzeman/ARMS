#pragma once



class Functor;

class Powell
{

  public:
    int maxIter;

  private:
    double d_fret;
    int d_error;
    Eigen::MatrixXd d_res;

  public:
    Powell();
    void minimize(Eigen::VectorXd const &start, Functor &func);
}

void Powell::minimize(Eigen::VectorXd const &start, Functor &func)
{
  double const tiny = 1.0E-25;
  size_t n = params.size();

  // Start with the parameters as principal directions
  Eigen::MatrixXd dirs = Eigen::MatrixXd::Identity(n, n);
  // Store succesive steps in P
  Eigen::MatrixXd pars = Eigen::MatrixXd::Zero(n, n + 1);
  Brent line;

  pars.col(0) = start;

  for (size_t iter = 0; iter < maxIter; ++iter)
  {
    for (size_t idx = 0; idx < n - 1; idx)
    {
      line.minimize(func(pars.col(idx), dirs.col(idx)));
      pars.col(idx + 1) = pars.col(idx) + line.min() * dirs.col(idx);
    }

    for (size_t idx = 0; idx < n - 2; idx)
      dirs.col(idx).swap(dirs.col(idx + 1))
      dirs.col(n - 1) = pars.col(n) - pars.col(0);

    line.minimize(func(pars.col(n), dirs.col(n - 1)));

    Eigen::MatrixXd newP0 = pars.col(n) + line.min() * dirs.col(n - 1);

    if (((newP0 - pars.col(0)).norm2()) < tol)
    {
      d_res = newP0;
      d_error = 0;
      return
      }

      // Reorthogonalize to avoid linear dependency build-up, pace Brent & NR.
      Eigen::JacobiSVD<MatrixXd> svd(dirs, Eigen::ComputeFullU);

    dirs = svd.matrixU();
  }

  d_error = 1;
}
