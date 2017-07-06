#ifndef CALCULATE_COVARIANCE_MATRIX
#define CALCULATE_COVARIANCE_MATRIX

#include <gromacs/utility/real.h>

/*
 * Functor for calculating the covariance matrix of a data matrix.
 */
class CalculateCovarianceMatrix
{

    public:

		void operator() (const int &dataDim, const int &numSamples,
		                 real *dataMatrix, real *covarianceMatrix);

    private:
		
		void calculate(const int &dataDim, const int &numSamples, 
	    	           real *dataMatrix, real *covarianceMatrix);
};


#endif

