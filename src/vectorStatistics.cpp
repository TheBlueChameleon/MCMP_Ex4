/* vectorStatistics.
 * Provides basic support for the std::vector<doube> class.
 * Implements average over the values: avg(vector) and
 * standard deviation: stdev(vector).
 */

// ========================================================================= //
// dependencies

#include <vector>
#include <cmath>
#include <algorithm>

#include "vectorStatistics.hpp"
#include "globals.hpp"

#include <iostream>

// ========================================================================= //
// direct statistical measures on a sample

double avg( std::vector<double> & data, unsigned int from, unsigned int till ) {
  double sum = 0;
  
  if (from > data.size()) {return NAN;}
  if (till > data.size()) {till = data.size();}
  
  for (auto i=from; i<till; i++) {
    sum += data[i];
  }
  
  return sum / (till - from);
}
// ------------------------------------------------------------------------- //
double variance ( std::vector<double> & data, unsigned int from, unsigned int till ) {
  /*  note: this returns sigma^2
   *                         ____________________________
   *                    _.  / sum(i=1..N) (x_i - x_med)²
   *            sigma =  | / ----------------------------
   *                     |/            N - 1
   */
  
  double med = avg(data, from, till),
         num = 0;
  
  if (from > data.size()) {return NAN;}
  if (till > data.size()) {till = data.size();}
  
  for (auto i=from; i<till; i++) {
    num += (data[i] - med) * (data[i] - med);
  }
  
  return num /= ((till - from) - 1);
}
// ......................................................................... //
double meanSquaredError ( std::vector<double> & data, unsigned int from, unsigned int till ) {
  // almost like variance, but with no bias compensation (divide by N, not by N-1)
  
  double med = avg(data, from, till),
         num = 0;
  
  if (from > data.size()) {return NAN;}
  if (till > data.size()) {till = data.size();}
  
  for (auto i=from; i<till; i++) {
    num += (data[i] - med) * (data[i] - med);
  }
  
//   std::cout << "MSE mean      : " << med << std::endl;
//   std::cout << "MSE sum       : " << num << std::endl;
//   std::cout << "MSE divides by: " <<        till - from  << std::endl;
//   std::cout << "MSE reval     : " << num / (till - from) << std::endl;
  
  return num / (till - from);
}
// ------------------------------------------------------------------------- //
double stdev( std::vector<double> & data, unsigned int from, unsigned int till ) {
  return std::sqrt(variance(data, from, till));
}
// ========================================================================= //
// autocorrelation time features

double autocorrelationTime(std::vector<double> & data, unsigned int discard) {
  double reVal = 0.5;
  
  for (auto rho : autocorrelationFunc(data, discard)) {reVal += rho;}
  
  return reVal;
}
// ......................................................................... //
double autocorrelationTimeFromCovFunc( std::vector<double> & rhoFunc ) {
  double reVal = 0.5;
  
  for (auto rho : rhoFunc) {reVal += rho;}
  
  return reVal;
}
// ......................................................................... //
std::vector<double> && autocorrelationFunc(std::vector<double> & data, unsigned int discard) {
  unsigned int N = data.size() - discard;
  
  std::vector<double>                 // std::vector CTor has form (count, value) where count is #elements and value is assigned to all elements of the vector.
    yPlus  (N, 0), yMinus(N, 0),
    denom  (N, 0);
  static std::vector<double>          // return a reference to this. ==> static lifetime.
    autoCov;
    autoCov.clear();
    autoCov.resize(N, 0);             // static means CTor is invoked only once -- manually prepare initial state.
  
  
  for   (unsigned int t=0; t<N  ; t++) {
    denom[t] = 1.0 / (N - t);
    
    for (unsigned int i=0; i<N-t; i++) {
      yMinus[t] += data[i  +discard];
      yPlus [t] += data[i+t+discard];
    }
    
    yMinus[t] *= denom[t];
    yPlus [t] *= denom[t];
  }
  
  
  for   (unsigned int t=0; t<N  ; t++) {
    for (unsigned int i=0; i<N-t; i++) {
      autoCov[t] += (data[i+discard] - yMinus[t]) * (data[i+t+discard] - yPlus[t]);
    }
    
    if (autoCov[t] < 0) {
        autoCov[t] = 0;
      break;
    }
    
    autoCov[t] *= denom[t];
  }
  
  double norm = autoCov[0];
  for (auto & rho : autoCov) {rho /= norm;}
  
  return std::move(autoCov);
}
// ------------------------------------------------------------------------- //
double stderror ( std::vector<double> & data, double tau, double sigma ) {
  if (std::isnan(tau  )) {tau   = autocorrelationTime(data);}
  if (std::isnan(sigma)) {sigma = stdev              (data);}
  
  return std::sqrt(2 * tau / data.size()) * sigma;
}
// ========================================================================= //
// error estimates on secondary quantities that are the variance of a sample

/* Block: Error propagation 
 * 
 * used secondary quantity has form of a variance:
 *    Q = variance(data) = <x²> - <x>²
 * using effective function in Markov chain points x_i
 *    g(x²_i, x_i) = 1 * x²_i   -   2 * (x_mean) * x_i
 *    where x_i is for data[tau..N]
 * with this, (error on f)² = sigma²_f
 *    sigma²_f = variance(g)
 * 
 * Note that this includes the unequilibrated part of the data.
 */
std::vector<double> && propagation_variance_g (std::vector<double> & data, unsigned int discard) {
  unsigned int N = data.size() - discard;
  static std::vector<double>          // return a reference to this. ==> static lifetime.
    reVal;
    reVal.clear();
    reVal.resize(N, 0);             // static means CTor is invoked only once -- manually prepare initial state.
  
  double mean = avg(data, discard);
  
  for (unsigned int i=0; i<N; i++) {
    reVal[i] = (data[i+discard] * data[i+discard]) - (2 * mean * data[i+discard]);
  }
  
  return std::move(reVal);
}
// ......................................................................... //
double propagation_variance_error(std::vector<double> & data, double scaling, unsigned int discard, double tau) {
  std::vector g = propagation_variance_g(data, discard);
  if ( std::isnan(tau) ) {tau = autocorrelationTime(g);}
  return std::sqrt( 2 * tau / g.size() ) * scaling * stdev(g);
}
// ------------------------------------------------------------------------- //
double blocking_variacne_error   (std::vector<double> & data, double scaling, unsigned int blocks) {
  if (blocks == 0) {return NAN;}
  
  unsigned int  N         = data.size(),
                blockSize = N / blocks;      
                  /* this wastes up to [blocks] Markov chain points due to
                   * the floating point part being truncated in the conversion
                   * to int.
                   * The effect still discards only a negligible subsection of 
                   * the total data set since usually [blocks << data.size()].
                   */
  std::vector<double> blockVal;
                      blockVal.resize(blocks, 0);
  
  for (unsigned int blockID = 1; blockID < blocks; blockID++) {
    blockVal[blockID] = 
      scaling * variance( data, 
                           blockID      * blockSize    , 
                          (blockID + 1) * blockSize - 1
                        );
  }
  
  return std::sqrt(meanSquaredError(blockVal, 1) / (blocks - 1));
  
}
// ------------------------------------------------------------------------- //
double bootstrap_variance_error  (std::vector<double> & data, double scaling, unsigned int discardedPointsCount, double tau, unsigned int M) {
  if ( std::isnan(tau) ) {tau = autocorrelationTime(data);}
  
  std::vector<std::vector<double>> pseudosamples(M);
  std::vector<            double > Q            (M);
  unsigned int                     N = data.size() - discardedPointsCount;
  int                              i = 0;
  
  if (!RNG_initialized) {
    std::cerr << "RNG not initialized!" << std::endl;
    return NAN;
  }
  
  for (auto & sample : pseudosamples) {
    sample.resize(N / tau);
    
    for (auto & point : sample) {
      point = data[gsl_rng_uniform_int(RNG, N) + discardedPointsCount];
    }
    
    Q[i] = scaling * variance(sample);
    i++;
  }
  
  return std::sqrt(meanSquaredError(Q));
  
}
