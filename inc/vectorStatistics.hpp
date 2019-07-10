#ifndef VECTORSTATISTICS_H
#define VECTORSTATISTICS_H

// ========================================================================= //
// dependencies

#include <climits>
#include <cmath>
#include <vector>

// ========================================================================= //
// proc

// ------------------------------------------------------------------------- //
// direct statistical measures on a sample

double avg              ( const std::vector<double> & data, const unsigned int from = 0, unsigned int till = INT_MAX );
double stdev            ( const std::vector<double> & data, const unsigned int from = 0, unsigned int till = INT_MAX );
double variance         ( const std::vector<double> & data, const unsigned int from = 0, unsigned int till = INT_MAX );
double meanSquaredError ( const std::vector<double> & data, const unsigned int from = 0, unsigned int till = INT_MAX );

// ------------------------------------------------------------------------- //
// autocorrelation time features

std::vector<double> autocorrelationFunc           (std::vector<double> & data, unsigned int discard = 0);
double              autocorrelationTime           (std::vector<double> & data, unsigned int discard = 0);
double              autocorrelationTimeFromCovFunc(std::vector<double> & rhoFunc);
double              stderror                      (std::vector<double> & data, double tau = NAN, double sigma = NAN);
  // computes the standard error given pre-computed autocorrelation time tau and standard deviation sigma.
  // if these two are NAN, the missing parameters are recomputed based on the content of data.

// ------------------------------------------------------------------------- //
// error estimates on secondary quantities that are the variance of a sample

std::vector<double> && propagation_variance_g    (std::vector<double> & data                    , unsigned int discard = 0);
double                 propagation_variance_error(std::vector<double> & data, double scaling = 1, unsigned int discard = 0, double tau  = NAN);

double                 blocking_variacne_error   (std::vector<double> & data, double scaling = 1, unsigned int blocks = 21);

double                 bootstrap_variance_error  (std::vector<double> & data, double scaling = 1, unsigned int discard = 0, double tau  = NAN, 
                                                  unsigned int M = 1000);

#endif// VECTORSTATISTICS_H
