/* Ising.cpp
 * defines the class Ising, which represents the spin grid and the relevant
 * methods
 */

// ========================================================================= //
// dependencies

#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "globals.hpp"
#include "Ising.hpp"

// ========================================================================= //
// macros

// translate 2D coordinates into 1D array index, incorporating periodic boundary conditions
#define ARRIDX(x, y) ( \
    (y >= L  ?  (y-L)  :  (y < 0  ?  (L+y)  :  y)) * (L) + \
    (x >= L  ?  (x-L)  :  (x < 0  ?  (L+x)  :  x))         \
  )

// get 2D coordinates from 1D array index
#define IDX_X(idx) ((idx) % (L))
#define IDX_Y(idx) ((idx) / (L))

// get automatic list of nearest neighbours
#define NEIGHBOURS(x, y) {\
    ARRIDX(x+1, y  ), \
    ARRIDX(x-1, y  ), \
    ARRIDX(x  , y+1), \
    ARRIDX(x  , y-1)  \
  }

// test whether x is between lo and hi
#define BETWEEN(x, lo, hi) ((x) >= (lo) ? ((x) <= (hi) ? true : false) : false)

// ========================================================================= //
// CTor, DTor

void Ising::init(double p) {
  // initialize with uniformly distributed spins, p(up) = p, p(down) = 1 - p.
  
  for (auto i=0; i<V; i++) {
    this->gridpoints.push_back(   gsl_rng_uniform(RNG) > p  ?  -1  :  +1   );
    
    this->neighbours.emplace_back();
    for (auto n : NEIGHBOURS(  IDX_X(i), IDX_Y(i)  )) {
      this->neighbours[i].push_back(n);
    }
  }
  
} 

// ----------------------------------------------------------------------- //

Ising::Ising  (IsingStart LS) {
  
  switch (LS) {
    case IsingStart::COLD :
      this->init(0.0);
      break;
      
    case IsingStart::HOT :
      this->init(0.5);
      break;
  }
}

// ....................................................................... //

Ising::Ising  (double p) {this->init(p);}

// ----------------------------------------------------------------------- //

Ising::~Ising () {
  
}

// ========================================================================= //
// lattice manipulation
  
int Ising::flipEnergy(const int i) const {
#if defined(ISING_DEBUG)
  if (!BETWEEN(i, 0, V)) {
    hDebug << "flipEnergy: attempting to flip invalid index " << i << std::endl;
    return -1;
  }
#endif
  
  int reVal = 0;
  for (auto j : this->neighbours[i]) {
    reVal += this->gridpoints[j];
  }
  
  return 2 * reVal * this->gridpoints[i];
}

// ----------------------------------------------------------------------- //

int Ising::flip(const int i) {
#if defined(ISING_DEBUG)
  if (!BETWEEN(i, 0, V)) {
    hDebug << "flipEnergy: attempting to flip invalid index " << i << std::endl;
  }
#endif
  
  return 2 * (this->gridpoints[i] *= -1);
}

// ========================================================================= //
// thermodynamic properties

int    Ising::energy()        const {
  /* note: H(grid) = -J * sum_{neighbours j of i} (s_i * s_j)
   * no contribution of B-field, and implicit J=1
   */
  
  double reVal = 0;
  
  for   (auto i=0; i<V; i++) {
    for (auto j : this->neighbours[i]) {
      reVal -= this->gridpoints[i] * this->gridpoints[j];
    }
  }
  
  return reVal / 2.0;
}
// ....................................................................... //
int    Ising::magnetization() const {
  int reVal = 0;
  
  for (auto i=0; i<V; i++) {
    reVal += this->gridpoints[i];
  }
  
  return reVal;
}
// ----------------------------------------------------------------------- //
double Ising::energyDensity()        const {return static_cast<double>(this->energy())        / this->gridpoints.size();}
// ....................................................................... //
double Ising::magnetizationDensity() const {return static_cast<double>(this->magnetization()) / this->gridpoints.size();}

// ========================================================================= //
// visualize

void Ising::show(std::ostream & target, std::string separator) const {
  target << "\t";
  
  for   (int x=0; x<L; x++) {
    target << x << "\t";
  }
  target << std::endl;
  
  for   (int y=0; y<L; y++) {
    target << y << "\t";
    for (int x=0; x<L; x++) {
      target << this->gridpoints[ ARRIDX(x, y) ] << separator;
    }
    target << std::endl;
  }
}
