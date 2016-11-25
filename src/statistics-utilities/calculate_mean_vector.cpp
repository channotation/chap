#include <iostream>
#include <algorithm>

#include <cblas.h>

#include "statistics-utilities/calculate_mean_vector.hpp"


/*
 * Operator overload for functor signature.
 */
void 
CalculateMeanVector::operator()(const int &dataDim,
                                const int &numSamples,
								real *dataMatrix,
								real *meanVector)
{
	this -> calculate(dataDim, numSamples, dataMatrix, meanVector);
}


/*
 * Function to calculate the mean of each row in a data matrix.
 */
void
CalculateMeanVector::calculate(const int &dataDim,
                               const int &numSamples,
         		               real *dataMatrix,
		                       real *meanVector)
{
    // fill mean vector with all zeros:
    std::fill_n(meanVector, dataDim, 0.0f);

    // sum over samples in each dimension:
    for(int i = 0; i < numSamples; i++)
    {
        cblas_saxpy(dataDim, 1.0f, &dataMatrix[i], numSamples, 
                    meanVector, 1);
    }

    // scale by 1/N to get maximum likelihood estimate of mean:
    cblas_sscal(dataDim, 1.0f/numSamples, meanVector, 1);
}
