/* Ising.hpp
 * defines the class Ising, which represents the spin grid and the relevant
 * methods
 */

#ifndef ISING_HPP
#define ISING_HPP

// ========================================================================= //
// dependencies

#include <vector>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "globals.hpp"

// ========================================================================= //
// symbols

enum class IsingStart {
  COLD, HOT
};

inline double IsingStartToP (IsingStart IS);

// ========================================================================= //
// class

class Ising {
private:
  /* TODO: transform globals L, V into properties of this class.
   * 
   */
  
  std::vector<int>                       gridpoints;
  std::vector<std::vector<unsigned int>> neighbours;
  
  double pStart = 0.0;
  double T      = 0.1;
  double lookupExp[4];
  
  std::vector<double> historyE;
  std::vector<double> historyM;
  
  double tauE   = NAN;
  double valE   = NAN;
  double errE   = NAN;
  double valC   = NAN;
  double errC   = NAN;
  
  double tauM   = NAN;
  double valM   = NAN;
  double errM   = NAN;
  double valX   = NAN;
  double errX   = NAN;
  
  void init (double p);
  
public:
  // ....................................................................... //
  // CTor, DTor
  
  Ising  (double T, IsingStart IS );
  Ising  (double T, double p = 0.0);
  ~Ising ();
  
  // ....................................................................... //
  // lattice manipulation
  
  int flipEnergy(const int i) const;
    // returns the energy difference due to a spin flip at grid index i 
    // performs boundary check on i, if the symbol LATTICE_DEBUG is defined
  
  int flip(const int i);
    // actually flips the spin at position i and returns the new spin state.
    // performs boundary check on i, if the symbol LATTICE_DEBUG is defined
  
  // ....................................................................... //
  // thermodynamic properties
    // xxxDensity refers to the property normalized to the grid volume
  
  int    energy()        const;
  int    magnetization() const;
  
  double energyDensity()        const;
  double magnetizationDensity() const;
  
  // ....................................................................... //
  // runtime parameters
  
  double getT() const;
  void   setT(double T);
  
  // ....................................................................... //
  // simulation
  
  void run  (IsingStart IS );
  void run  (double p = 0.0);
  
  void reset(IsingStart IS );
  void reset(double p = 0.0);
  
  // ....................................................................... //
  // read results
  
  const std::vector<double> & getHistoryE () const;
  const std::vector<double> & getHistoryM () const;
  
  double getTauE          ();
  double getTauM          ();
  
  double getValE          ();
  double getErrE          ();
  
  double getValM          ();
  double getErrM          ();
  
  double getValC          ();
  double getErrC          ();    // implements bootstrap
  
  double getValX          ();
  double getErrX          ();    // implements bootstrap
  
  
  // ....................................................................... //
  // output
  
  void dumpGrid    (std::ostream & target = std::cout, std::string separator = "\t") const;
  void dumpHistoryE(std::ostream & target = std::cout, std::string separator = "\t") const;
  void dumpHistoryM(std::ostream & target = std::cout, std::string separator = "\t") const;
  void dumpReport  (std::ostream & target = std::cout, std::string separator = "\t") const;
};

#endif//ISING_HPP
