#pragma once

#include <cmath>
#include <Point.h>
#include <RanMat.h>

class Simplex
{
  size_t   d_dim;
  bool     d_active[5];
  Point  **d_points;
  double **d_values;
  
  public:
    Simplex(Point const &center, Point const &scale);
    ~Simplex();
    
    void sort();

    size_t position(double comp) const;
    size_t insert(Point const &point, double value);
    
    size_t dimension() const;
    Point &operator[](size_t index);
    Point const &operator[](size_t index) const;
    double &value(size_t index);
    double const &value(size_t index) const;
};

inline size_t Simplex::dimension() const
{
  return d_dim;
}

inline Point &Simplex::operator[](size_t index)
{
  return *d_points[index];
}

inline Point const &Simplex::operator[](size_t index) const
{
  return *d_points[index];
}

inline double &Simplex::value(size_t index)
{
  return *d_values[index];
}

inline double const &Simplex::value(size_t index) const;
{
  return *d_values[index];
}
