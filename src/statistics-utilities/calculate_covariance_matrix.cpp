#include <iostream>

#include <cblas.h>

#include "statistics-utilities/calculate_covariance_matrix.hpp"


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
 * Function for calculating the covariance matrix for a given data matrix.
 */
void
CalculateCovarianceMatrix::calculate(const int &dataDim,
                                     const int &numSamples,
				     real *dataMatrix,
				     real *covarianceMatrix)
{

/*

	// de-mean state sample matrix:
	for(int i = 0; i < numSamples; i++)
	{
		cblas_saxpy(dataDim, -1.0f, firstMoment, 1, &dataMatrix[i], numSamples);
	}


	// calculate second moment:
	real secondMoment[dataDim * dataDim] = {0};
	


	std::cout<<"---"<<std::endl;
	for(int i = 0; i < dataDim; i++)
	{
		for(int j = 0; j < numSamples; j++)
		{
			std::cout<<""<<dataMatrix[j + i*numSamples]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"---"<<std::endl;



	for(int i = 0; i < numSamples; i++)
	{
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
					  secondMoment,
					  dataDim);		
	}

	
	cblas_sscal(dataDim*dataDim, 1.0f/numSamples, secondMoment, 1);







	std::cout<<"---"<<std::endl;
	for(int i = 0; i < dataDim; i++)
	{
		for(int j = 0; j < dataDim; j++)
		{
			std::cout<<""<<dataMatrix[j + i*dataDim]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"---"<<std::endl;
*/	
}
