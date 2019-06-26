 /* globals.hpp
 * behavioural parameters and global macros, debug file and GSL RNG
 */

// ========================================================================= //

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

// ========================================================================= //
// dependencies

#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// ========================================================================= //
// symbols

// ------------------------------------------------------------------------- //
// simulation behaviour

// system size
#define L 64

// system volume
#define V (L*L)

// steps in the MC simulation
#define N_MC  100000

// equilibration time
#define N_EQ  10000

#define REPORT_DIR "./txtout/"

// ------------------------------------------------------------------------- //
// meta behaviour
#define DEBUG_FILENAME "MCMP_2_Debug.txt"

#define DEBUG_LAUNCH hDebug << "Running Simulation. Launch time: " << __TIME__ << std::endl << SEPARATOR
#define DEBUG_END    hDebug << SEPARATOR << "Simulation ended in regular fashion. End time: " << __TIME__ << std::endl

#define DEBUG_FLUSH  hDebug << std::flush

#define ISING_DEBUG

// #define FEEDBACK_ONSCREEN

// ------------------------------------------------------------------------- //
// output shorthands

#define SMALLSEP  "// ......................................................................... //\n"
#define MEDSEP    "// ------------------------------------------------------------------------- //\n"
#define SEPARATOR "// ========================================================================= //\n"
#define PRINT_SMALLSEP   std::cout << SMALLSEP
#define PRINT_MEDSEP     std::cout << MEDSEP
#define PRINT_SEPARATOR  std::cout << SEPARATOR

#define FORMAT_INT    std::cout << std::showpos << std::fixed << std::setprecision(0)
#define FORMAT_DOUBLE std::cout << std::showpos << std::fixed << std::setprecision(5)
#define FORMAT_STD    std::cout << std::defaultfloat << std::setprecision(-1) << std::noshowpos

// ========================================================================= //
// global variables

extern std::ofstream    hDebug;
extern bool      RNG_initialized;
extern gsl_rng * RNG;

// ========================================================================= //
// procs

void RNG_init ();

#include <sstream>

// ========================================================================= //
// template procs

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

#endif//GLOBALS_HPP

