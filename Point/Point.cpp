#include <Point/Point.h>

inline Point &Point::operator=(Point const &other)
{
  sigma = other.sigma;
  m = other.m;
  a6 = other.a6;
  a7 = other.a7;
  a8 = other.a8;
  
  return *this;
}

inline Point &Point::operator+=(Point const &other)
{
  sigma += other.sigma;
  m     += other.m;
  a6    += other.a6;
  a7    += other.a7;
  a8    += other.a8;
  
  return *this;
}

inline Point &Point::operator-=(Point const &other)
{
  sigma -= other.sigma;
  m     -= other.m;
  a6    -= other.a6;
  a7    -= other.a7;
  a8    -= other.a8;
  
  return *this;
}

inline Point &Point::operator*=(double fac)
{
  sigma *= fac;
  m     *= fac;
  a6    *= fac;
  a7    *= fac;
  a8    *= fac;
  
  return *this;
}