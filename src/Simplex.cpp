#include <Simplex.h>
#include <iomanip>
#include <Log.h>

Simplex::Simplex(Data &data, Params &params, Weight w)
: d_dim(0), d_prec(params.prec), d_comp(data, params), d_values(0), d_weight(w)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  
  if (Log::ionode)
    log() << "Setting up initial simplex, using: " << std::endl;
  for (size_t idx = 0; idx < 6; ++idx)
    d_active[idx] = (std::abs(params.scale.coord[idx] / params.center.coord[idx]) > 1e-6);
  for (size_t idx = 0; idx < 6; ++idx)
    d_dim += d_active[idx] ? 1 : 0;
  
  d_points = new Point*[d_dim];
  std::fill_n(d_points, d_dim, static_cast< Point* >(0));
  d_values = 0;
  
  // Span the system from lower to upper for all active components
  d_points[0] = new Point(params.center);
  Point subt(params.scale);
  subt /= (d_dim - 1.0);
  *d_points[0] -= subt;
  size_t actIdx = 0;
  for (size_t idx = 1; idx < d_dim; ++idx)
  {
    while (!d_active[actIdx])
      ++actIdx;
    d_points[idx] = new Point(*d_points[0]);
    d_points[idx]->coord[actIdx] += params.scale.coord[actIdx];
    ++actIdx;
  }
  for (size_t idx = 0; idx < d_dim; ++idx)
    d_points[idx]->nonNegative();
  
  if (Log::ionode)
  {
    log() << *this << std::endl;
  // We want our values precise enough, so let's request things an order or magnitude
  // better than we would like this minimization to be.
    log() << "Calculating associated values.\n" << std::endl;
  }
  recalculate();  
  if (Log::ionode)
    log() << "Ordered simplex ready!\n" << std::endl;
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

void Simplex::recalculate()
{
  if (!d_values)
  {
    d_values = new double*[d_dim];
    for (size_t idx = 0; idx < d_dim; ++idx)
      d_values[idx] = new double;
  }
  for (size_t idx = 0; idx < d_dim; ++idx)
    *d_values[idx] = getVal(*d_points[idx]);

  sort();
  calcCenterOfGravity();
}

double Simplex::getVal(Point const &point)
{
  switch (d_weight)
  {
    case KOL:
      return d_comp.kolmogorov(point);
    case AVE:
      return d_comp.averages(point);
    default:
      return (-1.0);
  }
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
  while ((pos < d_dim) && (*d_values[pos] < val))
    ++pos;
  return pos;
}

size_t Simplex::constructProposal(double coeff)
{
  construct(coeff);
  if (Log::ionode)
    log() << "Proposing " << d_proposal << std::endl;
  d_propValue = getVal(d_proposal);
  return position(d_propValue);
}

bool Simplex::improveProposal(double coeff)
{
  Point oldProp = d_proposal;
  double oldVal = d_propValue;
  
  construct(coeff);
  if (Log::ionode)
    log() << "Proposing " << d_proposal << " as a further improvement." << std::endl;
  d_propValue = getVal(d_proposal);
  
  if (oldVal < d_propValue)
  {
    d_proposal = oldProp;
    d_propValue = oldVal;
    return false;
  }
  return true;
}

void Simplex::acceptProposal()
{
  *d_values[d_dim - 1] = d_propValue;
  *d_points[d_dim - 1] = d_proposal;
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

void Simplex::construct(double coeff)
{
  d_proposal  = d_cog;
  d_proposal -= *d_points[d_dim - 1];
  d_proposal *= coeff;
  d_proposal += d_cog;
  d_proposal.nonNegative();
}

void Simplex::reduceSimplex(double coeff)
{
  for (size_t idx = 1; idx < d_dim; ++idx)
  {
    *d_points[idx] -= *d_points[0];
    *d_points[idx] *= coeff;
    *d_points[idx] += *d_points[0];
    if (Log::ionode)
      log() << "Reduction of point " << idx << " would produce " << *d_points[idx] << std::endl;
    *d_values[idx] = getVal(*d_points[idx]);
  }
  sort();
}

std::ostream &operator<<(std::ostream &out, Simplex const &simplex)
{
  out << "        Sigma         m             a6            a7            a8            " << (simplex.d_values ? "|   D           " : "") << std::endl
      << "   ===========================================================================" << (simplex.d_values ? "================" : "") << std::endl;

  for (size_t idx = 0; idx < simplex.d_dim; ++idx)
  {
    out << "   " << idx << ".   " << std::setw(14) << std::left << simplex.d_points[idx]->coord[0]
                                  << std::setw(14) << std::left << simplex.d_points[idx]->coord[1]
                                  << std::setw(14) << std::left << simplex.d_points[idx]->coord[2]
                                  << std::setw(14) << std::left << simplex.d_points[idx]->coord[3]
                                  << std::setw(14) << std::left << simplex.d_points[idx]->coord[4];
    if (simplex.d_values)
      out << "|   "  << std::setw(10) << std::left << *simplex.d_values[idx];
    out << std::endl;
  }
  return out;
}

void Simplex::setWeight(Weight w)
{
  if (d_weight != w)
  {
    d_weight = w;
    recalculate();
  }
}
