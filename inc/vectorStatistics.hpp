#ifndef VECTORSTATISTICS_H
#define VECTORSTATISTICS_H

#include <climits>

double avg        ( std::vector<double> & data, unsigned int from = 0, unsigned int till = INT_MAX );
double stdev      ( std::vector<double> & data, unsigned int from = 0, unsigned int till = INT_MAX );
double variation  ( std::vector<double> & data, unsigned int from = 0, unsigned int till = INT_MAX );

#endif// VECTORSTATISTICS_H
