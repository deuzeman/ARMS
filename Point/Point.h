#pragma once

struct Point
{
  double coord[5];
 
  Point();
  Point(double sigma, double m, double a6, double a7, double a8);
  Point(Point const &other);
  
  Point &operator=(Point const &other);
  Point &operator+=(Point const &other);
  Point &operator-=(Point const &other);
  Point &operator*=(Point const &other);
}  

inline Point::Point(Point const &other)
{
  std::copy(other.coord, other.coord + 5, coord); 
}

Point &Point::operator=(Point const &other)
{
  if (&other != this)
    std::copy(other.coord, other.coord + 5, coord); 
  return *this;
}
