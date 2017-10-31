#include <gtest/gtest.h>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"


/*!
 * \brief Test fixture for the HistogramDensityEstimator.
 */
class AmiseOptimalBandwidthEstimatorTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        AmiseOptimalBandwidthEstimatorTest()
        {
            // parameters of normal distribution:
            real mu = -1.0;
            real sd = 3.0;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distribution(mu, sd);

            // create a random sample:
            size_t numSamples = 1e4;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
};


/*
 *
 */
TEST_F(AmiseOptimalBandwidthEstimatorTest, 
       AmiseOptimalBandwidthEstimatorGaussianDensityDerivativeFunctionalTest)
{
    // error tolerance for this test:
    real eps = std::numeric_limits<real>::epsilon();
 
    // pick an arbitrary standard deviation:
    real sigma = 3.3;

    // set the true values:
    real truePhi6
    real truePhi8

    // create bandwidth estimator:
    AmiseOptimalBandwidthEstimator bwe;

    // estimate the density derivative functionals:
    real phi6 = bwe.gaussianDensityDerivativeFunctional(sigma, 6);
    real phi8 = bwe.gaussianDensityDerivativeFunctional(sigma, 8);

    // assert correct value:
    ASSERT_NEAR(truePhi6, phi6, eps);
    ASSERT_NEAR(truePhi8, phi8, eps);
}
