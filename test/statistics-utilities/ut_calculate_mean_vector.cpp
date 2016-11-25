#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h> 

#include <gtest/gtest.h>

#include "statistics-utilities/calculate_mean_vector.hpp"


/*
 *
 */
class CalculateMeanVectorTest : public ::testing::Test
{
	protected:

		const static int seed_ = 15011991;

};


/*
 *
 */
TEST_F(CalculateMeanVectorTest, deterministicDataTest)
{
    // set parameters:
    int dataDim = 2;
    int numSamples = 3;

    // set up data matrix:
    real dataMatrix[dataDim * numSamples] = {-1.0f, 0.0f, 1.0f,
                                              1.0f, 2.0f, 3.0f};

    // calculate mean vector:
    real mv[dataDim];
    CalculateMeanVector calculateMeanVector;
    calculateMeanVector(dataDim, numSamples, dataMatrix, mv);
    
    // check that result is correct in each dimension:
    ASSERT_FLOAT_EQ(mv[0], 0.0f);
    ASSERT_FLOAT_EQ(mv[1], 2.0f);
}


/*
 *
 */
TEST_F(CalculateMeanVectorTest, randomDataTest)
{
    // set tolerance threshold:
    real relTol = 1e-2;

    // set parameters:
    int dataDim = 3;
    int numSamples = 1e5;
    
    real mean1 = -1.0f;
    real mean2 = 1.0f;
    real mean3 = 1e7f;
    
    real lo1 = -1.0f;
    real lo2 = -1.0f;
    real lo3 = -1.0f;
    
    real hi1 = 1.0f;
    real hi2 = 1.0f;
    real hi3 = 1.0f;

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

    // calculate mean vector:
    real meanVector[dataDim];
    CalculateMeanVector calculateMeanVector;
    calculateMeanVector(dataDim, numSamples, dataMatrix, meanVector);

    // check that mean is approximately correct for each dimension:
    ASSERT_NEAR((meanVector[0] - mean1)/mean1, 0.0, relTol);
    ASSERT_NEAR((meanVector[1] - mean2)/mean2, 0.0, relTol);
    ASSERT_NEAR((meanVector[2] - mean3)/mean3, 0.0, relTol);

    // clean up:
    delete[] dataMatrix;
}
