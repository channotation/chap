#include <iostream>

#include <cblas.h>

#include "statistics-utilities/calculate_covariance_matrix.hpp"
#include "statistics-utilities/calculate_mean_vector.hpp"


/*
 * Operator overload for functor signature.
 */
void CalculateCovarianceMatrix::operator()(const int &dataDim,
                                           const int &numSamples,
					   real *dataMatrix,
					   real *covarianceMatrix)
{
    this -> calculate(dataDim, numSamples, dataMatrix, covarianceMatrix);
}


/*
 * Function for calculating the covariance matrix for a given data matrix. Note
 * that consistent with Vanderbilt and Louie, this is a biased estimate of the 
 * covariance matrix, where each entry is scaled by 1/N rather than 1/(N-1).
 */
void
CalculateCovarianceMatrix::calculate(const int &dataDim,
                                     const int &numSamples,
								     real *dataMatrix,
								     real *covarianceMatrix)
{	
	// calculate the mean / first moment:
	real firstMoment[dataDim];
	CalculateMeanVector calculateMeanVector;
	calculateMeanVector(dataDim, numSamples, dataMatrix, firstMoment);
	
	// de-mean state sample matrix:
	for(int i = 0; i < numSamples; i++)
	{
		cblas_saxpy(dataDim, -1.0f, firstMoment, 1, &dataMatrix[i], numSamples);
	}
	
	// calculate second moment / covariance matrix:
	std::fill_n(covarianceMatrix, dataDim*dataDim, 0.0f);
	cblas_sgemm(CblasRowMajor,
				CblasNoTrans,
                CblasTrans,
				dataDim,
				dataDim,
				numSamples,
				1.0f,
				dataMatrix,
				numSamples,
				dataMatrix,
				numSamples,
				1.0f,
				covarianceMatrix,
				dataDim);		
	cblas_sscal(dataDim*dataDim, 1.0f/numSamples, covarianceMatrix, 1);
}
