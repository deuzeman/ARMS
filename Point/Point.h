#pragma once

struct Point
{
  double sigma;
  double m;
  double a6;
  double a7;
  double a8;
 
  Point();
  Point(double sigma, double m, double a6, double a7, double a8);
  Point(Point const &other);
  
  Point &operator=(Point const &other);
  Point &operator+=(Point const &other);
  Point &operator-=(Point const &other);
  Point &operator*=(Point const &other);
};

inline Point::Point()
 : sigma(0.0), m(0.0), a6(0.0), a7(0.0), a8(0.0)
{}

inline Point::Point(double sigma, double m, double a6, double a7, double a8)
  : sigma(sigma), m(m), a6(a6), a7(a7), a8(a8)
{}

inline Point::Point(Point const &other)
  : sigma(other.sigma), m(other.m), a6(other.a6), a7(other.a7), a8(other.a8)
{}
