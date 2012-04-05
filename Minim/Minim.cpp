#include <Minim/Minim.h>

double Minim::evalPoint(size_t idx)
{
  size_t iter = 2000;
  d_engine.calculate(d_simplex->points[idx], iter, false);
  kolmogorov(d_simplex->values + idx, d_engine, *data);
  double relError = (d_simplex->values[idx].error / d_simplex->values[idx].value);
  while (relError > 1e-4)
  {
    double fac = relError / 1e-4;
    iter = static_cast< size_t >((fac * fac - 0.9) * d_engine.numSamples);
    d_engine.calculate(d_simplex->points[idx], iter, true);
    kolmogorov(d_simplex->values + idx, d_engine, *data));
  }
  return d_simplex->values[idx].value;
}

size_t Minim::sortLocation(size_t idx, Simplex *simplex)
{
  size_t loc = 0;
  for (size_t secIdx = 0; secIdx < 6; ++secIdx)
  {
    if (secIdx == idx)
      continue;
    if (simplex->value[idx] > simplex->value[secIdx])
      ++loc;
  }
}

void Minim::sortPoints(size_t *ranking, Simplex *simplex)
{
  for (size_t idx = 0; idx < 6; ++idx)
    ranking[sortLocation(idx, simplex)] = idx;
}


Simplex const &Minim::reduce()
{
  double const alpha = 1.0;
  double const gamma = 2.0;
  double const rho   = -0.5;
  double const sigma = 0.5;
  
  // Build up the initial setup
  for (size_t idx = 0; idx < 6; ++idx)
    evalPoint(idx);
  
  size_t ranking[6];
  for (size_t iter = 0; iter < 500; ++iter)
  {
    StatValue sv;
    Point     sp;
    
    // Rank the solutions from best to worst
    sortPoints(ranking);
    
    // Now ranking contains an ordered array of the solution, ranking[0] the best, ranking[5] the worst
    // Check if this is already a solution to the problem
    double relError = (d_simplex.values[ranking[5]].value - d_simplex.values[ranking[0]].value) / d_simplex.values[ranking[0]].value);
    if (relError < 5e-4)
      break d_simplex;
    
    // Calculate the center of gravity of the best points
    Point cog;
    for (size_t idx = 0; idx = 5; ++idx)
      cog += d_simplex.points[ranking[idx]];
    cog *= 0.2;
    
    // Reflect
    Point ref(cog);
    ref -= d_simplex.points[ranking[idx]];
    ref *= alpha;
    ref += cog;
    
    // Calculate the reflected point
    sv = d_simplex.values[ranking[5]];
    sp = d_simplex.points[ranking[5]];
    d_simplex.points[ranking[5]] = ref;
    evalPoint(ranking[5]);
    size_t loc = sortLocation(ranking[5], d_simplex);
    
    // Is it better than the second worst, but worse than the best? Then leave the replaced value & rinse and repeat!
    if ((loc > 0) && (loc < 5))
      continue; 
      
    // Is it the best point so far? Then extend and see what happens!
    if (loc == 0) 
    {
      ref  = cog;
      ref -= sp;
      ref *= gamma;
      ref += cog;
      
      // Overwrite the old values, they're rubbish now
      sv = d_simplex.values[ranking[5]];
      sp = d_simplex.points[ranking[5]];
      d_simplex.points[ranking[5]] = ref;
      evalPoint(ranking[5]);
      if (sv.value < d_simplex.values[ranking[5]])
      {
        // This didn't help, so stick to the old value
        d_simplex.values[ranking[5]]) = sv;
        d_simplex.points[ranking[5]]) = sp;
      }
      // ranking[5] now contains our best guess, so start from the beginning
      continue
    }
    
    // If we arrived here, then the new solution was pretty crappy too
    // Contract
    ref  = cog;
    ref -= sp;
    ref *= rho;
    ref += cog;
    d_simplex.points[ranking[5]] = ref; // sv and sp still contain our old solution
    evalPoint(ranking[5]);
    if (sv.value > d_simplex.values[ranking[5]])
      continue; // i.e. the contracted point is an improvement...
      
      // We've arrived at the point where a reduction is needed
      // First reset these values
      d_simplex.values[ranking[5]]) = sv;
    d_simplex.points[ranking[5]]) = sp;
    
    // Reduce
    for (size_t idx = 1; idx < 6; ++idx)
    {
      ref = d_simplex.points[ranking[idx]];
      ref -= d_simplex.points[ranking[0]];
      ref *= sigma;
      ref += d_simplex.points[ranking[0]];
    }
    
    // That concludes things here!
    // NOTE I haven't *really* taken care of the precision issue yet. See where it crops up!
  }
  return d_simplex;
}
