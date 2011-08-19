#pragma once



class Brent
{
  public:
    double tol;
    double maxIter;
    double zeps;

    double ax;
    double bx;
    double cx;

  private:
    double d_fmin;
    double d_x;
    int d_error;

  public:
    Brent();
    void minimize(double (*func)(double const par));
	
	double const &min() const;
	double const &fmin() const;
};

inline Brent::Brent()
{}

inline double const &Brent::min() const
{
  return d_x;
}

inline double const &Brent::fmin() const
{
  return d_fmin;
}

void Brent::minimize(double (*func)(double const par))
{
  double const cgold = 0.3819660;

  double a = ((ax < cx) ? ax : cx);
  double c = ((cx > ax) ? cx : ax);
  double d = 0.0;
  double x = bx;
  double w = bx;
  double v = bx;

  double fx = func(x);
  double fv = fx;
  double fw = fx;

  d_error = -1;

  for (int iter = 0; iter < maxIter; ++iter)
  {
    double xm =   0.5 * (a + b);
    double tol1 = tol * std::abs(x) + zeps;
    double tol2 = 2.0 * tol1;

    if (std::abs(x - xm) <= (tol2 - 0.5 * (b -a)))
    {
      d_fmin = fx;
      d_xmin = x;
      d_error = 0;
	  return;
    }

    if (std::abs(e) > tol1)
    {
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);

      if (q > 0.0)
        p = -p;
      q = std::abs(q);
      double e_temp = e;
      e = d;
      if (abs(p) >= abs(0.5 * q * e_temp) || p <= q * (a - x) || p >= q * (b - x))
      {
        e = (x >= xm) ? a - x : b - x;
        d = cgold * x;
        d = cgold * e;
      }
      else
      {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = (xm >= x) ? tol : -tol;
      }
    }
    else
    {
      e = (x >= xm) ? a - x : b - x);
      d = cgold * e; 
    }

    u = (abs(d) >= tol1) ? x + d : x + ((d >= 0) ? tol1, -tol1);
    fu = func(u);

    if (fu <= fx)
    {
      if (u >= x)
        a = x;
      else
        b = x;
      v  =  w;
      w  =  x;
      x  =  u;
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
      if (fu <= fw || w == x)
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else
      {
        if (fu <= fv || v == x || v == w)
        {
          v = u;
          fv = fu;
        }
      }
    }
  }
  d_error = 1;
}
