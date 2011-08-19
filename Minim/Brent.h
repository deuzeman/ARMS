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
