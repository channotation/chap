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
       AmiseOptimalBandwidthEstimatorGaussianDensityRandomSampleTest)
{
    AmiseOptimalBandwidthEstimator bwe;

    // parameters of normal distribution:
    std::vector<real> mean = {-1.0, 0.0, 1.0, 1000.0, std::sqrt(2.0)};
    std::vector<real> sd = {1.0, 2.0, 3.0, 1.0, 0.5};
    std::vector<real> num = {10.0, 200.0, 1000.0, 50.0, 100};

    // prepare random distribution:
    std::default_random_engine generator;

    // generate test data sets:
    std::vector<std::vector<real>> testData;
    for(size_t i = 0; i < mean.size(); i++)
    {
        // prepare distribution:
        std::normal_distribution<real> distribution(mean[i], sd[i]);

        // draw random sample:
        std::vector<real> data;
        for(size_t j = 0; j < num[i]; j++)
        {
            data.push_back( distribution(generator) );
        }

        // estimate bandwidth:
        real bw = bwe.estimate(samples);

    }
}



















