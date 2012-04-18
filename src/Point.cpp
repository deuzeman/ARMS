#include <Point.h>

Point::Point()
{
  std::fill_n(coord, 5, 0.0);
}

Point::Point(double sigma, double m, double a6, double a7, double a8)
{
  coord[0] = sigma;
  coord[1] = m;
  coord[2] = a6;
  coord[3] = a7; 
  coord[4] = a8;
}

Point &Point::operator+=(Point const &other)
{
  for (size_t idx = 0; idx < 5; ++idx)
    coord[idx] += other.coord[idx];
  return *this;
}

Point &Point::operator-=(Point const &other)
{
  for (size_t idx = 0; idx < 5; ++idx)
    coord[idx] -= other.coord[idx];
  return *this;
}

Point &Point::operator*=(double fac)
{
  for (size_t idx = 0; idx < 5; ++idx)
    coord[idx] *= fac;
  return *this;
}

Point &Point::operator/=(double fac)
{
  for (size_t idx = 0; idx < 5; ++idx)
    coord[idx] /= fac;
  return *this;
}

Point &Point::operator=(Point const &other)
{
  if (&other != this)
    std::copy(other.coord, other.coord + 5, coord); 
  return *this;
}
