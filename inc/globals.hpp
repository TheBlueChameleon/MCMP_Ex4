 /* globals.hpp
 * behavioural parameters and global macros, debug file and GSL RNG
 */

// ========================================================================= //

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

// ========================================================================= //
// dependencies

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// ========================================================================= //
// DEFINE symbols

#define ISING_DEBUG
#define FEEDBACK_ONSCREEN

// ------------------------------------------------------------------------- //
// debug macros

#define DEBUG_LAUNCH hDebug << "Running Simulation. Launch time: " << __TIME__ << std::endl << SEPARATOR
#define DEBUG_END    hDebug << SEPARATOR << "Simulation ended in regular fashion. End time: " << __TIME__ << std::endl

#define DEBUG_FLUSH  hDebug << std::flush

// ------------------------------------------------------------------------- //
// output shorthands

#define SMALLSEP  "// ......................................................................... //\n"
#define MEDSEP    "// ------------------------------------------------------------------------- //\n"
#define SEPARATOR "// ========================================================================= //\n"

// ========================================================================= //
// global variables

// ------------------------------------------------------------------------- //
// simulation behaviour

extern const double outerT_lo;
extern const double outerT_hi;
extern const double outerT_dT;

extern const double innerT_lo;
extern const double innerT_hi;
extern const double innerT_dT;

// ------------------------------------------------------------------------- //
// log behaviour

extern const std::string   REPORT_DIR;
extern const std::string   DEBUG_FILENAME;
extern       std::ofstream hDebug;

// ------------------------------------------------------------------------- //
// GSL RNG

extern bool      RNG_initialized;
extern gsl_rng * RNG;

// ========================================================================= //
// procs

void RNG_init ();

// ========================================================================= //
// template procs

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

// ------------------------------------------------------------------------- //

// test whether x is between lo and hi
inline bool between(const double x, const double lo, const double hi) {return (x >= lo)  ?  ((x <= hi) ? true : false) : false;}

#endif//GLOBALS_HPP

