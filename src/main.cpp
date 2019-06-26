/* main.cpp
 * runs the Ising model simulation.
 */

// ========================================================================= //
 
#include <iostream>
#include <iomanip>
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
  
  // ----------------------------------------------------------------------- //
  // Metropolis-Hastings algorithm
  
  PRINT_SEPARATOR;
  std::cout << "Setup with dimension L=" << L << std::endl;
  std::cout << "\t#Monte Carlo sweeps              : " << N_MC << std::endl;
  std::cout << "\t#MC sweeps used for equilibration: " << N_EQ << std::endl;
  
  double lookupExp[3];
    /* we only need certain values exp(-i/T); however, by including exp(-0/T) 
     * in the table, we can use the return value of Latice::flipEnergy() as 
     * index of this C-array.
     */
  
  int    idx;               // index of the grid point to be flipped
  int    deltaE, totalE;    // change in energy due to the flip about to be done; total energy
  int            totalM;    // same for magnetization
  
  std::vector<double> historyE, historyM;
  std::ofstream       hReport;
  std::string         nReport;
  
  for   (auto T : {2.0, 2.3, 2.6}) {
    // ..................................................................... //
    // initialize system
    
    PRINT_SEPARATOR;
    std::cout << "initializing: T=" << T << std::endl;
    
    for (auto i=0; i<3; i++) {
      lookupExp[i] = std::exp(-4.0*i/T);      // TODO: idx -> exp(??)
    }
    
    for (auto pStart : {0.0, 0.5, 0.75}) {
      PRINT_MEDSEP;
      std::cout << "initial probability for spin-up: pStart=" << pStart << std::endl;
      
      Ising grid = Ising(pStart);
      
      totalE = grid.energy();
      totalM = grid.magnetization();
      
      historyE.clear();
      historyM.clear();
      
      historyE.push_back(         static_cast<double>(totalE) / V) ;
      historyM.push_back(std::abs(static_cast<double>(totalM) / V));
      
      nReport  = REPORT_DIR;
      nReport += "Report_T=";
      nReport += to_string_with_precision(     T, 1);
      nReport += "_pStart=";
      nReport += to_string_with_precision(pStart, 1);
      nReport += ".txt";
      
      
      // ..................................................................... //
      // run MCMH
      
      for   (int i=0; i<N_MC; i++) {
        for (int sweep=0; sweep<V; sweep++) {
          idx = gsl_rng_uniform_int(RNG, V);
          deltaE = grid.flipEnergy(idx);
          
          if (deltaE <= 0) {
            totalM += grid.flip(idx);
            totalE += deltaE;
            
          } else {
            double randVal = gsl_rng_uniform(RNG);
            
            if (randVal < lookupExp[deltaE/4]) {
              totalM += grid.flip(idx);
              totalE += deltaE;
            }
          }
          
        }
        
        historyE.push_back(         static_cast<double>(totalE) / V) ;
        historyM.push_back(std::abs(static_cast<double>(totalM) / V));
        
        if (! (i% 1000)) {
          std::cout << "." << std::flush;
        }
      }
      
      std::cout << "done." << std::endl;
      
      std::cout << "\theat density                   : " <<             avg      (historyE, N_EQ)     << std::endl;
      std::cout << "\tmagnetization density          : " <<             avg      (historyM, N_EQ)     << std::endl;
      std::cout << "\tSpecific heat density          : " << (1/(T*T)) * variation(historyE, N_EQ) * V << std::endl; 
      std::cout << "\tMagnetic susceptibility density: " << (1/  T  ) * variation(historyM, N_EQ) * V << std::endl;
      
      
      // ..................................................................... //
      // output to file
      
      hReport.open(nReport);
      hReport << std::fixed << std::setprecision(5);
      hReport << "# ID, energy density, magnetization density" << std::endl;
      
      for (auto i=0; i<N_MC; i++) {
//          if (! (i% 10)) {        // reduce the number of plotted data points by factor %__
          hReport << i << "\t" << historyE[i] << "\t" << historyM[i] << std::endl;
//          }
      }
      
      hReport.close();
    }
  }
  
  // ----------------------------------------------------------------------- //
  // tidy up
  
  std::cout << "done." << std::endl;
  DEBUG_END;
}
