#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>

struct Point
{
  double coord[5];
 
  friend std::ostream &operator<<(std::ostream &out, Point const &point);
  
  Point();
  Point(double sigma, double m, double a6, double a7, double a8);
  Point(Point const &other);
  
  Point &operator=(Point const &other);
  Point &operator+=(Point const &other);
  Point &operator-=(Point const &other);
  Point &operator*=(double fac);
  Point &operator/=(double fac);
  
  void nonNegative();
};

inline Point::Point(Point const &other)
{
  std::copy(other.coord, other.coord + 5, coord); 
}

inline std::ostream &operator<<(std::ostream &out, Point const &point)
{
  out << '[' << point.coord[0] << ' ' << point.coord[1] << ' ' << point.coord[2] << ' ' << point.coord[3] << ' ' << point.coord[4] << ']';
  return out;
}

inline void Point::nonNegative()
{
  for (size_t idx = 0; idx < 5; ++idx)
    if (coord[idx] < 0.0)
      coord[idx] = 0.0;
}
