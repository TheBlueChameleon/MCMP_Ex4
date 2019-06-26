/* Ising.hpp
 * defines the class Ising, which represents the spin grid and the relevant
 * methods
 */

#ifndef ISING_HPP
#define ISING_HPP

// ========================================================================= //
// dependencies

#include <vector>

#include <iostream>
#include <fstream>
#include <string>

#include "globals.hpp"

// ========================================================================= //
// symbols

enum class IsingStart {
  COLD, HOT
};

// ========================================================================= //
// class

class Ising {
private:
  std::vector<int>                       gridpoints;
  std::vector<std::vector<unsigned int>> neighbours;
  
  void init (double p);
  
public:
  // ....................................................................... //
  // CTor, DTor
  
  Ising  (IsingStart IS);
  Ising  (double p = 0.5);
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
  // Visualize
  
  void show (std::ostream & target = std::cout, std::string separator = "\t") const;
};

#endif//ISING_HPP
