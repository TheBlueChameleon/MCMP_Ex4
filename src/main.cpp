/* main.cpp
 * runs the Ising model simulation.
 */

// ========================================================================= //
 
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "globals.hpp"
#include "Ising.hpp"
#include "vectorStatistics.hpp"

// ========================================================================= //
// proc

int main () {
  // ----------------------------------------------------------------------- //
  // initialize
  
  DEBUG_LAUNCH;
  RNG_init();
  
  const std::string filename = "./out/evolution.dat";
  std::ofstream hf_evolution(filename);
  
  hf_evolution
      << "# "
      << "temperature T\t"
      << "heat density e\terror on e\t" 
      << "absolute magnetization density |m|\terror on m\t"
      << "specific heat density c\terror on c\t"
      << "magnetic susceptibility chi\terror on chi"
      << std::endl;
  
  // ----------------------------------------------------------------------- //
  // Metropolis-Hastings algorithm
  
  std::cout << SEPARATOR;
  std::cout << "Setup with dimension L=" << L << std::endl;
  std::cout << "\t#Monte Carlo sweeps              : " << N_MC << std::endl;
  
  std::ofstream       hReport;
  std::string         nReport;
  
  double dT = outerT_dT;
  Ising model(0.0, IsingStart::COLD);
  
  for (auto T = outerT_lo; T <= outerT_hi; T += dT) {
    if (between(T, innerT_lo, innerT_hi)) {dT = innerT_dT;}
    else                                  {dT = outerT_dT;}
    
    model.setT(T);
    model.run();
    
    // file out
    hf_evolution << T << "\t";
    
    if (std::isnan(model.getTauE()) || 
        std::isnan(model.getTauM())
    ) {
      hf_evolution << "# --- insufficient data to compute meaningfull data ---" << std::endl;
      continue;
    }
    
    hf_evolution << model.getValE() << "\t";
    hf_evolution << model.getErrE() << "\t";
    
    hf_evolution << model.getValM() << "\t";
    hf_evolution << model.getErrM() << "\t";
    
    hf_evolution << model.getValC() << "\t";
    hf_evolution << model.getErrC() << "\t";
    
    hf_evolution << model.getValX() << "\t";
    hf_evolution << model.getErrX();
    hf_evolution << std::endl;
    
    std::cout << "post files" << std::endl;
  }
  
  // ----------------------------------------------------------------------- //
  // tidy up
  
  hf_evolution.close();
  
  std::cout << "done." << std::endl;
  DEBUG_END;
}
