#ifndef MISC_STATS_UTILITIES_HPP
#define MISC_STATS_UTILITIES_HPP

#include <gromacs/utility/real.h>

/*
 * Function to calculate mean of floating point array.
 */
real
calculateArrayMean(const int &numSamples, real *dataArray);

#endif

