/* main.cpp
 * runs the Ising model simulation.
 */

// ========================================================================= //
// dependencies

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "globals.hpp"
#include "Ising.hpp"
#include "vectorStatistics.hpp"

// ========================================================================= //
// switches

// #define RUN_METROPOLIS
#define RUN_WOLFF

// ========================================================================= //
// proc

int main () {
  // ----------------------------------------------------------------------- //
  // initialize
  
  DEBUG_LAUNCH;
  RNG_init();
  
  const std::string fn_Metropolis = "./out/results_Metropolis.dat";
  const std::string fn_Wolff      = "./out/results_Wolff.dat";
  std::ofstream     fh_Metropolis(fn_Metropolis);
  std::ofstream     fh_Wolff     (fn_Wolff);
  
  const unsigned int L = 64;
  
#ifdef RUN_METROPOLIS
  fh_Metropolis
      << "# "
      << "temperature T\t"
      << "heat density e\terror on e\t" 
      << "absolute magnetization density |m|\terror on m\t"
      << "specific heat density c\terror on c\t"
      << "magnetic susceptibility chi\terror on chi\t"
      << "average cluster size"
      << std::endl;
#endif
  
#ifdef RUN_WOLFF
  fh_Wolff
      << "# "
      << "temperature T\t"
      << "heat density e\terror on e\t" 
      << "absolute magnetization density |m|\terror on m\t"
      << "specific heat density c\terror on c\t"
      << "magnetic susceptibility chi\terror on chi"
      << std::endl;
#endif
  
  // ----------------------------------------------------------------------- //
  // Metropolis-Hastings algorithm
  
  std::cout << SEPARATOR;
  std::cout << "Setup with dimension L=" << L << std::endl;
  
  double dT = outerT_dT;
  Ising model(L, IsingStart::COLD);
  
  for (auto T = outerT_lo; T <= outerT_hi; T += dT) {
    if (between(T, innerT_lo, innerT_hi)) {dT = innerT_dT;}
    else                                  {dT = outerT_dT;}
    
#ifdef RUN_METROPOLIS
    model.setT(T);
    model.run_Metropolis();
    
    std::cout << "computing primary and secondary quantities and their errors...";
    
    // file out
    fh_Metropolis << T << "\t";
    
    if (std::isnan(model.getTauE()) || 
        std::isnan(model.getTauM())
    ) {
      fh_Metropolis << "# --- insufficient data to compute meaningfull data ---" << std::endl;
      goto skipPointMetropolis;
    }
    
    fh_Metropolis << model.getValE() << "\t";
    fh_Metropolis << model.getErrE() << "\t";
    
    fh_Metropolis << model.getValM() << "\t";
    fh_Metropolis << model.getErrM() << "\t";
    
    fh_Metropolis << model.getValC() << "\t";
    fh_Metropolis << model.getErrC() << "\t";
    
    fh_Metropolis << model.getValX() << "\t";
    fh_Metropolis << model.getErrX();
    fh_Metropolis << std::endl;
    
    
    skipPointMetropolis:
    std::cout << "done" << std::endl;
#endif
    
    
#ifdef RUN_WOLFF
    model.setT(T);
    model.run_Wolff();
    
    std::cout << "computing primary and secondary quantities and their errors..." << std::flush;
    
    // file out
    fh_Wolff << T << "\t";
    
    if (std::isnan(model.getTauE()) || 
        std::isnan(model.getTauM())
    ) {
      fh_Wolff << "# --- insufficient data to compute meaningfull data ---" << std::endl;
      goto skipPointWolff;
    }
    
    std::cout << std::endl;
    std::cout << "tau_E: " << model.getTauE() << std::endl;
    std::cout << "tau_M: " << model.getTauM() << std::endl;
    
    fh_Wolff << model.getValE() << "\t";
    fh_Wolff << model.getErrE() << "\t";
    
    fh_Wolff << model.getValM() << "\t";
    fh_Wolff << model.getErrM() << "\t";
    
    fh_Wolff << model.getValC() << "\t";
    fh_Wolff << model.getErrC() << "\t";
    
    fh_Wolff << model.getValX() << "\t";
    fh_Wolff << model.getErrX();
    fh_Wolff << std::endl;
    
    skipPointWolff:
    std::cout << "done" << std::endl;
#endif
  }
  
  // ----------------------------------------------------------------------- //
  // tidy up
  
  fh_Metropolis.close();
  fh_Wolff     .close();
  
  std::cout << "simulation done." << std::endl;
  DEBUG_END;
}
