#ifndef CALCULATE_MEAN_VECTOR
#define CALCULATE_MEAN_VECTOR

#include <gromacs/utility/real.h>

/*
 * Functor for calculating the vector of means for a given data matrix.
 */
class CalculateMeanVector
{
	public:

		void operator() (const int &dataDim, const int &numSamples,
		                 real *dataMatrix, real *meanVector);

	private:

		void calculate(const int &dataDim, const int &numSamples,
		               real *dataMatrix, real *meanVector);

};


#endif

