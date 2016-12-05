#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h> 

#include <gtest/gtest.h>

#include <lapacke.h>

#include "statistics-utilities/calculate_covariance_matrix.hpp" 


/*
 * Test fixture for testing the covariance matrix calculation functor. 
 */
class CalculateCovarianceMatrixTest : public ::testing::Test
{
	protected:

		const static int seed_ = 15011993;

};


/*
 * Tests covariance matrix calculation on small, low-dimensional data set of 
 * known values and makes sure all covariance matrix elements agree with the 
 * correct values (of a biased covariance matrix). There are both a correlated 
 * and an anti-correlated case.
 */
TEST_F(CalculateCovarianceMatrixTest, deterministicDataTest)
{
	// set up parameters:
	int dataDim = 2;
	int numSamples = 3;

	// set up data structure and functor:
	real covarianceMatrix[dataDim * dataDim];
	CalculateCovarianceMatrix calculateCovarianceMatrix;

	// set up data matrix for perfectly correlated case:
	real dataMatrix1[dataDim * numSamples] = {-1.0f, 0.0f, 1.0f,
	                                          -1.0f, 0.0f, 1.0f};
	
	// calculate covariance matrix:	
	calculateCovarianceMatrix(dataDim, numSamples, dataMatrix1, covarianceMatrix);

	// check that result is correct:
	ASSERT_FLOAT_EQ(covarianceMatrix[0], 2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[1], 2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[2], 2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[3], 2.0f/3.0f);	
	
	// set up data matrix for perfectly anti-correlated case:
	real dataMatrix2[dataDim * numSamples] = {-1.0f, 0.0f,  1.0f,
	                                           1.0f, 0.0f, -1.0f};
	
	// calculate covariance matrix:	
	calculateCovarianceMatrix(dataDim, numSamples, dataMatrix2, covarianceMatrix);

	// check that result is correct:
	ASSERT_FLOAT_EQ(covarianceMatrix[0], 2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[1], -2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[2], -2.0f/3.0f);
	ASSERT_FLOAT_EQ(covarianceMatrix[3], 2.0f/3.0f);	
}


/*
 * Tests the calculation of the covariance matrix on a medium sized, low
 * dimensional data set of random numbers. It is asserted that the covariance
 * is symmatric and that it is positive semi-definite, i.e. that its eigen-
 * values are all larger than or equal to zero.
 */
TEST_F(CalculateCovarianceMatrixTest, randomDataTest)
{
    // set parameters:
    int dataDim = 3;
    int numSamples = 1e2;
    
    real mean1 = -1.0f;
    real mean2 = 1.4f;
    real mean3 = 1.e5f;
    
    real lo1 = -1.0f;
    real lo2 = -1e0f;
    real lo3 = -1.0f;
    
    real hi1 = 1.0f;
    real hi2 = 1e4f;
    real hi3 = 1.0f;

	// set up data structures and functor:
	real covarianceMatrix[dataDim * dataDim];
	CalculateCovarianceMatrix calculateCovarianceMatrix;
	
	// set up random number generator:
    gmx::DefaultRandomEngine rng(seed_);                     
    gmx::UniformRealDistribution<real> distribution1(lo1, hi1);
    gmx::UniformRealDistribution<real> distribution2(lo2, hi2);  
    gmx::UniformRealDistribution<real> distribution3(lo3, hi3); 

    // set up data matrix:
    real *dataMatrix = new real[dataDim * numSamples];
    for(int i = 0; i < numSamples; i++)
    {
		dataMatrix[i + 0*numSamples] = distribution1(rng) + mean1;
		dataMatrix[i + 1*numSamples] = distribution2(rng) + mean2;
		dataMatrix[i + 2*numSamples] = distribution3(rng) + mean3;

    }

	// calculate covariance matrix:	
	calculateCovarianceMatrix(dataDim, numSamples, dataMatrix, covarianceMatrix);

	// check symmetry of matrix:
	for(int i  = 0; i < dataDim; i++)
	{
		for(int j = 0; j < dataDim; j++)
		{
			ASSERT_FLOAT_EQ(covarianceMatrix[i*dataDim + j], 
			                covarianceMatrix[j*dataDim + i]);
		}
	}


	// calculate eigenvalues of from upper triangle of covariance matrix:
	real ev1[dataDim] = {0};
	int status1 = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 
	              	            'N', 
							    'U', 
							    dataDim,
							    covarianceMatrix,
							    dataDim, 
							    ev1);

	// check that eigenvalue calculation terminated successfully:
	ASSERT_EQ(0, status1);

	// calculate eigenvalues of from lower triangle of covariance matrix:
	real ev2[dataDim] = {0};
	int status2 = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 
	              	            'N', 
							    'L', 
							    dataDim,
							    covarianceMatrix,
							    dataDim, 
							    ev2);

	// check that eigenvalue calculation terminated successfully:
	ASSERT_EQ(0, status2);

	// check eigenvalues:	
	for(int i = 0; i < dataDim; i++)
	{
		// covariance matrix should be positive-semidefinite:	
		ASSERT_GE(ev1[i], 0.0);
		ASSERT_GE(ev2[i], 0.0);
	}
}

