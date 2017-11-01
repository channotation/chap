#include <random>

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


/*!
 *
 */
TEST_F(AmiseOptimalBandwidthEstimatorTest, 
       AmiseOptimalBandwidthEstimatorGaussianDensityNonnegativityTest)
{
/*
    // parameters of normal distribution:
    std::vector<real> mean = {-1.0, 0.0, 1.0, 1000.0, std::sqrt(2.0)};
    std::vector<real> sd = {1.0, std::sqrt(2.0), 100.0, 1.0, 2.0};
    std::vector<int> num = {100, 10, 10, 10, 10};

    // prepare random distribution:

    // generate test data sets:
    for(size_t i = 0; i < mean.size(); i++)
    {
        // prepare distribution:
        std::default_random_engine generator;
        std::normal_distribution<real> distribution(mean[i], sd[i]);

        std::cout<<"mu = "<<mean[i]<<"  "
                 <<"sd = "<<sd[i]<<"  "
                 <<"N = "<<num[i]<<"  ";

        // draw random sample:
        std::vector<real> data;
        for(size_t j = 0; j < num[i]; j++)
        {
            data.push_back( distribution(generator) );
            std::cout<<data.back()<<"  ";
        }

        // estimate bandwidth:
        AmiseOptimalBandwidthEstimator bwe;
        real bw = bwe.estimate(data);

        std::cout<<"bw = "<<bw<<std::endl;

        // make sure bandwidth is positive:
        ASSERT_LT(0.0, bw);
    }
   */ 
}


/*!
 *
 */
TEST_F(AmiseOptimalBandwidthEstimatorTest, 
       AmiseOptimalBandwidthEstimatorGaussianDensityReferenceImplementationTest)
{
    // error tolerance:
    real eps = std::numeric_limits<real>::epsilon();

    // randomly generated sample data (from standard normal):
    std::vector<real> sample = {-0.4462099,
                                -1.5673473,
                                -0.4568042,
                                -0.4100268,
                                 0.7060560,
                                -0.5724193,
                                 0.8842685,
                                -1.8202730,
                                 0.6218274,
                                -0.5386946};

    std::vector<real> test;
    for(int i = 0; i < 1000; i++)
    {
        for(auto s : sample)
        {
            test.push_back(s);
        }
    }

    // estimate AMISE-optimal bandwidth for this:
    AmiseOptimalBandwidthEstimator bwe;
    real bw = bwe.estimate(test);

    // compare to reference implementation value of bandwidth:
    real trueBw = 0.156558;
    ASSERT_NEAR(trueBw, bw, eps);
}















