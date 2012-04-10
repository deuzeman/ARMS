#include <Simplex.h>

Simplex::Simplex(Point const &center, Point const &scale)
: d_dim(0)
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
  
  for (size_t idx = 0; idx < d_dim; ++idx)
  {
    d_values[idx] = new double;
    *d_values[idx] = d_ranmat.kolmogorov(*d_points[idx]);
  }
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
}

size_t Simplex::position(double val) const
{
  size_t pos = 0;
  while ((*d_values[pos] < val) && (val < d_dim))
    ++pos;
  return pos;
}

size_t Simplex::insert(Point const &point)
{
  double value = d_ranmat.kolmogorov(point);
  size_t pos = position(value);
  if (pos <= *d_values[d_dim - 1])
  {
    *d_values[d_dim - 1] = value;
    *d_points[d_dim - 1] = point;
    sort();
  }
  return pos;
}
