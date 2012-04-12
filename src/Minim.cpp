#include <Minim.h>
#include <Log.h>

Simplex const &Minim::reduce()
{
  log() << "\nStarting minimization routine!\n\n";
  for (size_t iter = 0; iter < 500; ++iter)
  {
    log() << "At iteration " << iter << ":\n" << d_simplex;
    // Check if this is already a solution to the problem
    if (d_simplex.converged())
    {
      log() << "Simplex values indicate convergence!\nFinal solution:\n" << d_simplex;
      break;
    }
    
    Point proposed;
    
    log() << "REFLECT\n";
    size_t loc = d_simplex.constructProposal(d_alpha); // Reflect
    
    // Is it the best point so far? Then extend and see what happens!
    // Either way, we've made a major step and we can drop the worst value.
    // We can just drop through to the next if-case
    if (loc == 0) 
    {
      log() << "EXTEND\n";
      d_simplex.improveProposal(d_gamma); // Extend
    }
    
    // Is our proposal better than the second worst? That'll do us fine, get rid of the worst one and repeat
    if (loc < (d_simplex.dimension() - 1))
    {
      log() << "Reflection accepted.\n";
      d_simplex.acceptProposal();
      continue; 
    }
    else
      log() << "Reflection rejected.\n";

    
    // If we arrived here, then the proposal was pretty crappy too
    // That is, it was *perhaps* better than the worst know point.
    // At this point, we contract and don't care about the old proposal either way.
    
    log() << "CONTRACT\n";
    loc = d_simplex.constructProposal(d_rho); // Contract

    if (loc > d_simplex.dimension())
    {
      log() << "Contraction accepted.\n";
      d_simplex.acceptProposal();
      continue; // i.e. the contracted point is an improvement, so we continue on
    }
    else
      log() << "Contraction rejected.\n";
      
    // We've arrived at the point where a reduction is needed
    // Reduce
    log() << "REDUCE\n";
    d_simplex.reduceSimplex(d_sigma);
    
    // That concludes things here!
    // NOTE I haven't *really* taken care of the precision issue yet. See where it crops up!
  }
  return d_simplex;
}
