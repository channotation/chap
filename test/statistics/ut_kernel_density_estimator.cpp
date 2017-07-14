#include <algorithm>
#include <limits>
#include <random>

#include <gtest/gtest.h>

#include "statistics/kernel_density_estimator.hpp"


/*
 *
 */
class KernelDensityEstimatorTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        KernelDensityEstimatorTest()
        {
            // parameters of normal distribution:
            real mu = -1.0;
            real sd = 3.0;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distribution(mu, sd);

            // create a random sample:
            size_t numSamples = 1;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
};


/*!
 * This test checks that the evaluation points are calculated correctly.
 * Specifically they are required to be well ordered and to cover the complete
 * data range. In addition, it is checked that the number of evaluation points
 * is a power of two. All checks are performed for a range of different
 * maximum evaluation point distances.
 */
TEST_F(KernelDensityEstimatorTest, KernelDensityEstimatorEvalPointTest)
{
    // create kernel density estimator and set parameters:
    KernelDensityEstimator kde;
    kde.setBandWidth(1.0);
    kde.setEvalRangeCutoff(0.0);

    // define range of evaluation point distances:
    std::vector<real> maxEvalPointDistances = {10.0, 1.0, 0.1, 0.01, 1e-5};

    // perform tests for each value of evaluation distance:
    for(auto maxEvalPointDist : maxEvalPointDistances)
    {
        // set max eval dist to current value:
        kde.setMaxEvalPointDist(maxEvalPointDist);

        // obtain evaluation points:
        std::vector<real> evalPoints = kde.createEvaluationPoints(testData_);

        // assert that number of eval points is power of two:
        ASSERT_TRUE( (evalPoints.size() & (evalPoints.size() - 1)) == 0 );

        // assert that evaluation points are well ordered:
        for(size_t i = 0; i < evalPoints.size() - 1; i++)
        {
            ASSERT_LT(evalPoints.at(i), evalPoints.at(i+1));
        }

        // get actual data range:
        real dataLo = *std::min_element(testData_.begin(), testData_.end());
        real dataHi = *std::max_element(testData_.begin(), testData_.end());

        // assert that data range is covered:
        ASSERT_GE(dataLo, evalPoints.front());
        ASSERT_LE(dataHi, evalPoints.back());
    }
}


/*
 *
 */
TEST_F(KernelDensityEstimatorTest, KernelDensityEstimatorGaussianDensityTest)
{    
    // integral over density should be within eps of one:
    real eps = std::sqrt(std::numeric_limits<real>::epsilon());

    // define a number of bandwidths to test:    
    std::vector<real> bandWidths = {10.0, 1.0, 0.1, 0.01, 1e-5};

    // perform tests for all the above bandwidths:
    for(auto bw : bandWidths)
    {
        // create kernel density estimator and set parameters:
        // (note that large cutoff makes probability mass outside cutoff range
        // negligible).
        KernelDensityEstimator kde;
        kde.setBandWidth(bw);
        kde.setEvalRangeCutoff(10.0); 
        kde.setMaxEvalPointDist(0.1*bw); // if >> bw, integral is not one
        kde.setKernelFunction(eKernelFunctionGaussian);

        // get evaluation points and calculate density:
        std::vector<real> evalPoints = kde.createEvaluationPoints(testData_);
        std::vector<real> density = kde.calculateDensity(testData_, evalPoints);

        // assert non-negativity of density:
        for(auto d : density)
        {
            ASSERT_LE(0.0, d);
        }

        // get actual evaluation point distance:
        real evalPointDist = (evalPoints.back() - evalPoints.front());
        evalPointDist /= (evalPoints.size() - 1);

        // assert that density integrates to one:
        real integral = 0.0;
        for(auto d : density)
        {
            integral += d;
        }
        integral *= evalPointDist;
        ASSERT_NEAR(1.0, integral, eps);
    }
}

