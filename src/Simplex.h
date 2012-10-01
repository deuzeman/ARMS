#pragma once

#include <cmath>
#include <Comparator.h>
#include <Point.h>
#include <RanMat.h>

class Simplex
{
  size_t     d_dim;
  bool       d_active[5];
  Point    **d_points;
  double   **d_values;
  
  double     d_prec_a;
  double     d_prec_k;
  Comparator d_comp;
    
  Point      d_cog;
  Point      d_proposal;
  double     d_propValue;
  
  int        d_rank;
  
  friend std::ostream &operator<<(std::ostream &out, Simplex const &simplex);
  
  public:
    Simplex(Data &data, Params &params, int type = Comparator::AVE);
    ~Simplex();
    
    size_t constructProposal(double coeff);
    bool improveProposal(double coeff);
    void acceptProposal();
    
    void reduceSimplex(double coeff);
    
    size_t dimension() const;
    Point &operator[](size_t index);
    Point const &operator[](size_t index) const;
    double &value(size_t index);
    double const &value(size_t index) const;
    
    bool converged() const;
    
    void writeProposal() const;

    void setType(int type);
    int type() const;

  private:
    void construct(double coeff);
    size_t position(double comp) const;
    void sort();
    void calcCenterOfGravity();
    void recalculate();
    double getVal(Point const &point);
};

std::ostream &operator<<(std::ostream &out, Simplex const &simplex);

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

inline double const &Simplex::value(size_t index) const
{
  return *d_values[index];
}

inline void Simplex::writeProposal() const
{
  log() << d_proposal << std::endl;
}

inline double Simplex::getVal(Point const &point)
{
  return d_comp.deviation(point);
}

inline int Simplex::type() const
{
  return d_comp.type();
}


