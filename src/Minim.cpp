#include <Minim.h>
#include <Log.h>

Simplex const &Minim::reduce()
{
  if (Log::ionode)
    log() << "\nStarting minimization routine!\n\n";
  for (size_t iter = 0; iter < 500; ++iter)
  {
    if (Log::ionode) 
      log() << "At iteration " << iter << ":\n" << d_simplex;
    // Check if this is already a solution to the problem
    if (d_simplex.converged())
    {
      if (d_simplex.getWeight() == AVE)
      {
        if (Log::ionode)
          log() << "Simplex values indicate heuristic convergence using average value.\nSwitching to KS.\n";
        d_simplex.setWeight(KOL);
        if (Log::ionode)
          log() << "Under changed weighting, simplex has assumed the following form:\n" << d_simplex << std::endl;
        continue;
      }
      if (Log::ionode)
        log() << "Simplex values indicate convergence!\nFinal solution:\n" << d_simplex << std::endl;
      break;
    }
    
    Point proposed;
    if (Log::ionode)
      log() << "== REFLECT ==\n";
    size_t loc = d_simplex.constructProposal(d_alpha); // Reflect
    
    // Is it the best point so far? Then extend and see what happens!
    // Either way, we've made a major step and we can drop the worst value.
    // We can just drop through to the next if-case
    if (loc == 0) 
    {
      if (Log::ionode)
        log() << "== EXTEND ==\n";
      d_simplex.improveProposal(d_gamma); // Extend
    }
    
    // Is our proposal better than the second worst? That'll do us fine, get rid of the worst one and repeat
    if (loc < (d_simplex.dimension() - 1))
    {
      if (Log::ionode)
        log() << "Reflection accepted.\n";
      d_simplex.acceptProposal();
      continue; 
    }
    else
    {
      if (Log::ionode)
        log() << "Reflection rejected.\n";
    }
    
    // If we arrived here, then the proposal was pretty crappy too
    // That is, it was *perhaps* better than the worst know point.
    // At this point, we contract and don't care about the old proposal either way.
    if (Log::ionode)
      log() << "== CONTRACT ==\n";
    loc = d_simplex.constructProposal(d_rho); // Contract

    if (loc < d_simplex.dimension())
    {
      if (Log::ionode)
        log() << "Contraction accepted.\n";
      d_simplex.acceptProposal();
      continue; // i.e. the contracted point is an improvement, so we continue on
    }
    if (Log::ionode)
      log() << "Contraction rejected.\n";
      
    // We've arrived at the point where a reduction is needed
    // Reduce
    if (Log::ionode)
      log() << "== REDUCE ==\n";
    d_simplex.reduceSimplex(d_sigma);
    
    // That concludes things here!
    // NOTE I haven't *really* taken care of the precision issue yet. See where it crops up!
  }
  return d_simplex;
}
