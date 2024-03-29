/* Ising.cpp
 * defines the class Ising, which represents the spin grid and the relevant
 * methods
 */

// ========================================================================= //
// dependencies

#include <vector>
#include <array>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "globals.hpp"
#include "vectorStatistics.hpp"
#include "Ising.hpp"

// ========================================================================= //
// enum translator

inline double IsingStartToP (IsingStart IS) {
  switch (IS) {
    case IsingStart::COLD : return 0.0;
    case IsingStart::HOT  : return 0.5;
  }
  
  // this should not be able to happen
  std::cerr << "Warning: Encountered an invalid IsingStart symbol" << std::endl;
  return NAN;
}

// ========================================================================= //
// CTor, DTor

void Ising::init(unsigned int L, double p) {
  // initialize with uniformly distributed spins, p(up) = p, p(down) = 1 - p.
  
  this->pStart = p;
  this->L      = L;
  this->V      = L * L;
  
  this->gridpoints.resize(this->V);
  this->neighbours.resize(this->V);
  
  for (auto i=0u; i<V; i++) {
    auto theNeighbours = this->neighboursOf(arrayCoordinateX(i), arrayCoordinateY(i));
    std::copy(theNeighbours.cbegin(), theNeighbours.cend(), std::back_inserter(this->neighbours[i]));
  }
  
  this->reset();      // sets initial spins according to pStart
} 

// ----------------------------------------------------------------------- //

Ising::Ising  (unsigned int L, IsingStart IS) {this->init(L, IsingStartToP(IS));}
// ....................................................................... //
Ising::Ising  (unsigned int L, double pStart) {this->init(L, pStart           );}

// ----------------------------------------------------------------------- //
Ising::~Ising () {
  
}

// ========================================================================= //
// lattice geometry

// translate 1D array index into 2D coordinates
inline unsigned int Ising::arrayCoordinateX (const unsigned int idx) const {return idx % this->L;}
inline unsigned int Ising::arrayCoordinateY (const unsigned int idx) const {return idx / this->L;}

// translate 2D coordinates into 1D array index, incorporating periodic boundary conditions
inline unsigned int Ising::arrayIndex(const int x, const int y) const {
  return 
    (y >= static_cast<int>(this->L)  ?  y-this->L  :  (y < 0  ?  (this->L+y)  :  y)) * this->L + 
    (x >= static_cast<int>(this->L)  ?  x-this->L  :  (x < 0  ?  (this->L+x)  :  x));
}

// get automatic list of nearest neighbours
inline std::array<unsigned int, 4> Ising::neighboursOf (const int x, const int y) const {
  std::array<unsigned int, 4> reVal;
  
  reVal[0] = this->arrayIndex(x+1, y  );
  reVal[1] = this->arrayIndex(x  , y+1);
  reVal[2] = this->arrayIndex(x-1, y  );
  reVal[3] = this->arrayIndex(x  , y-1);
  
  return reVal;
}


// ========================================================================= //
// lattice manipulation
  
int Ising::flipEnergy(const int i) const {
#if defined(ISING_DEBUG)
  if (!between(i, 0, this->V)) {
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
  if (!between(i, 0, this->V)) {
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
  
  for   (unsigned int i=0; i<this->V; i++) {
    reVal -= this->gridpoints[i] * this->gridpoints[this->neighbours[i][0]];
    reVal -= this->gridpoints[i] * this->gridpoints[this->neighbours[i][1]];
  }
  
  return reVal;
}
// ....................................................................... //
int    Ising::magnetization() const {
  int reVal = 0;
  
  for (unsigned int i=0; i<this->V; i++) {
    reVal += this->gridpoints[i];
  }
  
  return reVal;
}
// ----------------------------------------------------------------------- //
double Ising::energyDensity()        const {return static_cast<double>(this->energy())        / this->gridpoints.size();}
// ....................................................................... //
double Ising::magnetizationDensity() const {return static_cast<double>(this->magnetization()) / this->gridpoints.size();}

// ========================================================================= //
// runtime parameters

double Ising::getT() const {return this->T;}
// ....................................................................... //
void   Ising::setT(const double T, const double pStart) {
  this->reset(pStart);
  
  this->T = T;
  
  for (auto i = 0u; i<4u; i++) {
    lookupExp[i] = std::exp(-4.0*i/T);
  }
}

// ========================================================================= //
// simulation

void Ising::run_Metropolis () {
  int    idx;               // index of the grid point to be flipped
  int    deltaE, totalE;    // change in energy due to the flip about to be done; total energy
  int            totalM;    // same for magnetization
  
#ifdef FEEDBACK_ONSCREEN
  std::cout << MEDSEP;
  std::cout << "running new Metropolis chain at T=" << this->T << " with p_up = " << pStart << std::endl;
#endif
  
  this->reset();
  
  totalE = this->energy();
  totalM = this->magnetization();
  
  this->historyE.push_back(         static_cast<double>(totalE) / this->V) ;
  this->historyM.push_back(std::abs(static_cast<double>(totalM) / this->V));
  
  
  for   (auto i     = 0u;     i < N_MC   ;     i++) {
    for (auto sweep = 0u; sweep < this->V; sweep++) {
      idx = gsl_rng_uniform_int(RNG, this->V);
      deltaE = this->flipEnergy(idx);
      
      if (deltaE <= 0) {
        totalM += this->flip(idx);
        totalE += deltaE;
        
      } else {
        double randVal = gsl_rng_uniform(RNG);
        
        if (randVal < lookupExp[deltaE/4]) {
          totalM += this->flip(idx);
          totalE += deltaE;
        }
      }
      
    }
    
    historyE.push_back(         static_cast<double>(totalE) / this->V) ;
    historyM.push_back(std::abs(static_cast<double>(totalM) / this->V));
    
    
#ifdef FEEDBACK_ONSCREEN
    if (! (i% 1000)) {
      std::cout << "." << std::flush;
    }
#endif
  }
  
#ifdef FEEDBACK_ONSCREEN
  std::cout << "done" << std::endl;
#endif
}
// ----------------------------------------------------------------------- //
void Ising::run_Wolff() {
  unsigned int idx;
  unsigned int cluster_ID = 0;
  int          seedspin;
  double       probability_to_add = 1 - std::exp(-2.0 / this->T);
  
  std::vector<unsigned int> cluster;
  std::vector<double      > lengths;
  
#ifdef FEEDBACK_ONSCREEN
  std::cout << MEDSEP;
  std::cout << "running new Wolff chain at T=" << this->T << " with p_up = " << pStart << std::endl;
#endif
  
  this->reset(pStart);
  
  this->historyE.push_back(this->energyDensity()) ;
  this->historyM.push_back(this->magnetizationDensity());
  
  for   (auto i = 0u; i < N_MC; i++) {
    cluster.clear();
    
    idx = gsl_rng_uniform_int(RNG, this->V);
    
    seedspin = this->gridpoints[idx];
    this->gridpoints[idx] *= -1;
    
    cluster.push_back(idx);
    cluster_ID = 0;
    
    do {
      for (auto neighbour : this->neighbours[cluster[cluster_ID]]) {
        if (this->gridpoints[neighbour] == seedspin) {
          double randVal = gsl_rng_uniform(RNG);
          
          if (randVal < probability_to_add) {
            this->gridpoints[neighbour] *= -1;
            cluster.push_back(neighbour);
          } 
        }
        
      }
      
      cluster_ID++;
    } while (cluster_ID < cluster.size());
    
    lengths.push_back(static_cast<double>(cluster.size()) / this->V);
    
    this->historyE.push_back(         this->energyDensity()        ) ;
    this->historyM.push_back(std::abs(this->magnetizationDensity()));
    
    
#ifdef FEEDBACK_ONSCREEN
    if (! (i% 1000)) {
      std::cout << "." << std::flush;
    }
#endif
  }
  
  this->cLen = avg(lengths);
  this->cErr = stderror(lengths);
  
#ifdef FEEDBACK_ONSCREEN
  std::cout << "done" << std::endl;
#endif
}
// ----------------------------------------------------------------------- //
void Ising::reset(IsingStart IS) {this->reset(IsingStartToP(IS));}
// ....................................................................... //
void Ising::reset(double p) {
  if (!std::isnan(p)) {this->pStart = p;}
  
  for (auto i=0u; i<this->V; i++) {
    this->gridpoints[i] = (gsl_rng_uniform(RNG) < pStart  ?  +1  :  -1);
  }
  
  this->historyE.clear();
  this->historyM.clear();
  
  this->valE = NAN;
  this->valC = NAN;
  this->tauE = NAN;
  this->errE = NAN;
  this->errC = NAN;
  
  this->valM = NAN;
  this->valX = NAN;
  this->tauM = NAN;
  this->errM = NAN;
  this->errX = NAN;
}

// ========================================================================= //
// read results

const std::vector<double> & Ising::getHistoryE() const {return this->historyE;}
// ....................................................................... //
const std::vector<double> & Ising::getHistoryM() const {return this->historyM;}
// ----------------------------------------------------------------------- //
double Ising::getTauE () {
  if (std::isnan(this->tauE)) {
    this->tauE = autocorrelationTime(this->historyE);
    if (std::isnan(this->tauE)) {return NAN;}
    
    this->tauE = autocorrelationTime(this->historyE, 20 * this->tauE);
    if (20 * this->tauE > this->historyE.size()) {this->tauE = NAN;}
  }
  
  return this->tauE;
}
// ....................................................................... //
double Ising::getTauM () {
  if (std::isnan(this->tauM)) {
    this->tauM = autocorrelationTime(this->historyM);
    if (std::isnan(this->tauM)) {return NAN;}
    
    this->tauM = autocorrelationTime(this->historyM, 20 * this->tauM);
    if (20 * this->tauM > this->historyM.size()) {this->tauM = NAN;}
  }
  
  return this->tauM;
}
// ----------------------------------------------------------------------- //
double Ising::getValE () {
  if (std::isnan(this->valE)) {
    this->valE = avg(this->historyE, 20 * this->getTauE());
  }
  return this->valE;
}
// ....................................................................... //
double Ising::getErrE () {
  if (std::isnan(this->errE)) {this->errE = stderror(historyE, this->getTauE());}
  return this->errE;
}
// ----------------------------------------------------------------------- //
double Ising::getValM () {
  if (std::isnan(this->valM)) {
    this->valM = std::abs(avg(this->historyM, 20 * this->getTauM()));
  }
  return this->valM;
}
// ....................................................................... //
double Ising::getErrM () {
  if (std::isnan(this->errM)) {this->errM = stderror(historyM, this->getTauM());}
  return this->errM;
}
// ----------------------------------------------------------------------- //
double Ising::getValC () {
  if (std::isnan(this->valC)) {
    this->valC = (1/(this->T*this->T)) * variance(this->historyE, 20 * this->getTauE()) * V;
  }
  return this->valC;
}
// ....................................................................... //
double Ising::getErrC () {
  if (std::isnan(this->errC)) {
    this->errC = bootstrap_variance_error  (this->historyE, V / (T*T), 20 * this->tauE, this->tauE);
  }
  return this->errC;
}
// ----------------------------------------------------------------------- //
double Ising::getValX () {
  if (std::isnan(this->valX)) {
    this->valX = (1/(    this->T    )) * variance(this->historyM, 20 * this->getTauM()) * V;
  }
  return this->valX;
}
// ....................................................................... //
double Ising::getErrX () {
  if (std::isnan(this->errX)) {
    this->errX = bootstrap_variance_error  (this->historyM, V / T, 20 * this->tauM, this->tauM);
  }
  return this->errX;
}
// ----------------------------------------------------------------------- //
double Ising::getCLen () const {return this->cLen;}
// ....................................................................... //
double Ising::getCErr () const {return this->cErr;}
// ========================================================================= //
// visualize

void Ising::dumpGrid(std::ostream & target, std::string separator) const {
  target << "\t";
  
  for   (unsigned int x=0; x<L; x++) {
    target << x << "\t";
  }
  target << std::endl;
  
  for   (unsigned int y=0; y<L; y++) {
    target << y << "\t";
    for (unsigned int x=0; x<L; x++) {
      target << this->gridpoints[ arrayIndex(x, y) ] << separator;
    }
    target << std::endl;
  }
}
