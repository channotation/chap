#include <random>

#include <gtest/gtest.h>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"

/*!
 * \brief Test fixture for the AmiseOptimalBandwidthEstimator.
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

        };
};


/*!
 * Checks the validity of the AMISE optimal bandwidth estimates by comparison
 * to the bandwidth estimated from Silverman's rule of thumb. This is carried
 * out for a Gaussian distribution (for which Silverman's rule should be exact
 * in the limit of large N) for a variety of sample sizes and parameters in the
 * original distribution.
 */
TEST_F(AmiseOptimalBandwidthEstimatorTest, 
       AmiseOptimalBandwidthEstimatorGaussianDensitySilvermanTest)
{
    // tolerance threshold:
    real tol = 0.1;

    // parameters:
    int numReps = 2;
    std::vector<real> mean = {-10, 0.0, std::sqrt(2.0)};
    std::vector<real> numSamples = {250, 500, 750};
    std::vector<real> standardDeviation = {100.0, 1.0, 0.01};

    // build sample from Guassian with varying standard deviation:
    for(auto mu : mean)
    {
        for(auto num : numSamples)
        {
            for(auto sd : standardDeviation)
            {
                // prepare Gaussian distribution:
                std::default_random_engine generator;
                std::normal_distribution<real> distribution(mu, sd);
                
                for(size_t i = 0; i < numReps; i++)
                {            
                    // draw a sample:
                    std::vector<real> sample;
                    for(int j = 0; j < num; j++)
                    {
                        sample.push_back( distribution(generator) );
                    }

                    // obtain sample standard deviation:
                    SummaryStatistics sumStat;
                    for(auto s : sample)
                    {
                        sumStat.update(s);
                    }
                    real ssd = sumStat.sd();
                
                    // estimate bandwidth:
                    AmiseOptimalBandwidthEstimator bwe;
                    real bw = bwe.estimate(sample);

                    // Silverman's rule of thumb should be accurate for pure Gaussian:
                    real silverman = 1.06*ssd/std::pow(num, 1.0/5.0);

                    // valid bandwidth should always be positive:
                    ASSERT_GT(bw, 0.0);
                    ASSERT_GT(silverman, 0.0);

                    // check consistency between Silverman and AMISE estimates:
                    ASSERT_NEAR(std::fabs(bw - silverman)/silverman, 0.0, tol);
                }
            }
        }
    }
}

