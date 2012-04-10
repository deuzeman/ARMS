#include <Simplex.h>

Simplex::Simplex(Point const &center, Point const &scale, Comparator &comp, double prec)
: d_dim(0), d_prec(prec), d_comp(comp)
{
  for (size_t idx = 0; idx < 5; ++idx)
    d_active[idx] = (std::abs(scale.coord[idx] / center.coord[idx]) > 1e-6);
  for (size_t idx = 0; idx < 5; ++idx)
    d_dim += d_active[idx] ? 1 : 0;
  
  d_points = new Point*[d_dim];
  d_values = new double*[d_dim];
  
  // Span the system from lower to upper for all active components
  d_points[0] = new Point(center);
  *d_points[0] -= scale;
  size_t actIdx = 0;
  for (size_t idx = 1; idx < d_dim; ++idx)
  {
    while (!d_active[actIdx])
      ++actIdx;
    d_points[idx] = new Point(*d_points[0]);
    d_points[idx]->d_coord[actIdx] += 2.0 * scale.d_coord[actIdx];
    ++actIdx;
  }
  
  // We want our values precise enough, so let's request things an order or magnitude
  // better than we would like this minimization to be.
  d_comp.setPrecision(0.1 * comp);
  for (size_t idx = 0; idx < d_dim; ++idx)
  {
    d_values[idx] = new double;
    *d_values[idx] = d_comp.kolmogorov(*d_points[idx]);
  }
  
  calcCenterOfGravity();
}

Simplex::~Simplex()
{
  for (size_t idx = 0; idx < d_dim; ++idx)
  {
    delete d_points[idx];
    delete d_values[idx];
  }
  delete d_points;
  delete d_values;
}

void Simplex::sort()
{ 
  for (size_t idx = 0; idx < d_dim - 1; ++idx)
  {
    size_t best = idx;
    double bestVal = *d_values[idx];
    
    for (size_t run = idx + 1; run < d_dim; ++run)
      if (*d_values[run] < bestVal)
      {
        best = run;
        bestVal = *d_values[run];
      }
    
    if (best != idx)
    {
      Point  *tp = d_points[idx];
      double *tv = d_values[idx];
      d_points[idx] = d_points[best];
      d_values[idx] = d_values[best];
      d_points[best] = tp;
      d_values[best] = tv;
    }
  }
  
  calcCenterOfGravity();
}

size_t Simplex::position(double val) const
{
  size_t pos = 0;
  while ((*d_values[pos] < val) && (val < d_dim))
    ++pos;
  return pos;
}

size_t Simplex::constructProposal(double coeff)
{
  construct(coeff);
  d_propValue = d_comp.kolmogorov(d_proposed);
  return position(d_propValue);
}

bool Simplex::improveProposal(double coeff)
{
  Point oldProp = d_proposal;
  double oldVal = d_propValue;
  
  construct(coeff);
  d_propValue = d_comp.kolmogorov(d_proposed);
  
  if (d_oldVal < d_propValue)
  {
    d_proposal = oldProp;
    d_propValue = oldval;
    return false;
  }
  return true;
}

void Simplex::acceptProposal()
{
  *d_values[d_dim - 1] = d_propValue;
  *d_points[d_dim - 1] = d_proposed;
  sort();
  calcCenterOfGravity(); // Might be often redundant, but it's cheap
}

void Simplex::calcCenterOfGravity()
{
  std::fill_n(d_cog.coord, 5, 0.0);
  for (size_t idx = 0; idx < d_dim - 1; ++idx)
    d_cog += *d_points[idx];
  d_cog *= (1.0 / (d_dim - 1));
}

void Simplex::construct(double coeff) const
{
  d_proposed  = d_cog;
  d_proposed -= d_simplex[d_dim - 1];
  d_proposed *= coeff;
  d_proposed += d_cog;
}

void Simplex::reduceSimplex(double coeff)
{
  for (size_t idx = 1; idx < d_dim; ++idx)
  {
    *d_points[idx] -= *d_points[0];
    *d_points[idx] *= d_sigma;
    *d_points[idx] += *d_points[0];
    *d_values[idx] = d_comp.kolmogorov(*d_points[idx]);
  }
}
