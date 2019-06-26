/* globals.hpp
 * behavioural parameters and global macros, debug file and GSL RNG
 */

// ========================================================================= //
// dependencies

#include <fstream>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "globals.hpp"

// ========================================================================= //
// debug

std::ofstream hDebug(DEBUG_FILENAME);

// ========================================================================= //
// GSL RNG

bool      RNG_initialized = false;
gsl_rng * RNG;

void RNG_quit ();

// ------------------------------------------------------------------------- //

void RNG_init () {
  gsl_rng_env_setup();
  RNG = gsl_rng_alloc (gsl_rng_default);    // default is MT19937
  
  time_t seed = std::time(nullptr);
  gsl_rng_set(RNG, seed);
  hDebug << "RNG of type " << gsl_rng_name (RNG) << " initialized with seed " << seed << std::endl;
  
  std::atexit(RNG_quit);
  
  RNG_initialized = true;
}

// ------------------------------------------------------------------------- //

void RNG_quit () {
  gsl_rng_free (RNG);
  RNG_initialized = false;
}
