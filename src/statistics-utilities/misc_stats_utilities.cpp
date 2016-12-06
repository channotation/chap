#include "statistics-utilities/misc_stats_utilities.hpp"

/*
 * Calculates mean of all values in an array.
 */
real
calculateArrayMean(const int &numSamples, real *dataArray)
{
	real mean = 0.0f;
	for(int i=0; i<numSamples; i++)
	{
		mean += dataArray[i];
	}
	mean /= numSamples;

	return(mean);
}

