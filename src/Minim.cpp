#include <Minim/Minim.h>

Simplex const &Minim::reduce()
{  
  size_t ranking[6];
  for (size_t iter = 0; iter < 500; ++iter)
  {
    // Now ranking contains an ordered array of the solution, ranking[0] the best, ranking[5] the worst
    // Check if this is already a solution to the problem
    if (d_simplex.converged())
      break;
    
    Point proposed;
    
    size_t loc = d_simplex.constructProposal(d_alpha); // Reflect
    
    // Is it the best point so far? Then extend and see what happens!
    // Either way, we've made a major step and we can drop the worst value.
    // We can just drop through to the next if-case
    if (loc == 0) 
      d_simplex.improveProposal(d_gamma); // Extend
    
    
    // Is our proposal better than the second worst? That'll do us fine, get rid of the worst one and repeat
    if (loc < (d_simplex.dimension() - 1))
    {
      d_simplex.acceptProposal();
      continue; 
    }
    
    // If we arrived here, then the proposal was pretty crappy too
    // That is, it was *perhaps* better than the worst know point.
    // At this point, we contract and don't care about the old proposal either way.
    
    loc = d_simplex.constructProposal(d_rho); // Contract
    if (loc > d_dim);
    {
      d_simplex.acceptProposal();
      continue; // i.e. the contracted point is an improvement, so we continue on
    }
      
    // We've arrived at the point where a reduction is needed
    // Reduce
    d_simplex.reduceSimplex(d_sigma);
    
    // That concludes things here!
    // NOTE I haven't *really* taken care of the precision issue yet. See where it crops up!
  }
  return d_simplex;
}
