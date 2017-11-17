#include <random>

#include <gtest/gtest.h>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"

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
            size_t numSamples = 1e0;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
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
    int numReps = 3;
    std::vector<real> mean = {-10, 0.0, std::sqrt(2.0)};
    std::vector<real> numSamples = {500, 1000, 1500};
    std::vector<real> standardDeviation = {100.0, 10.0, 1.0, 0.1, 0.01};

    // build sample from Guassian with varying standard deviation:
//    std::vector<real> bandWidth;
//    std::vector<real> sampleStandardDeviation;
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
    //                sampleStandardDeviation.push_back(ssd);

                    // shift and scale data:
                    GaussianDensityDerivative gdd;
                    auto ss = gdd.getShiftAndScaleParams(sample, sample);
                    gdd.shiftAndScale(sample, ss.first, ss.second);
                
                    // estimate bandwidth:
                    AmiseOptimalBandwidthEstimator bwe;
                    real bw = bwe.estimate(sample) / ss.second;
    //                bandWidth.push_back(bw);

                    // Silverman's rule of thumb should be accurate for pure Gaussian:
                    real silverman = 1.06*ssd/std::pow(num, 1.0/5.0);

                    // valid bandwidth should always be positive:
                    ASSERT_GT(bw, 0.0);
                    ASSERT_GT(silverman, 0.0);

    /*
                    std::cout.precision(5);
                    std::cout<<"n = "<<num<<"  "
                             <<"sd = "<<ssd<<"  "
                             <<"bw = "<<bw<<"  "
                             <<"silverman = "<<silverman<<"  "
                             <<"abs = "<<std::fabs(bw - silverman)<<"  "
                             <<"rel = "<<std::fabs(bw - silverman)/silverman<<"  "
                             <<"ratio = "<<silverman/bw<<"  "
                             <<std::endl;
    */
                    // check consistency between Silverman and AMISE estimates:
                    ASSERT_NEAR(std::fabs(bw - silverman)/silverman, 0.0, tol);
                }
            }
        }
    }


/*
    std::cout.precision(5);
    for(size_t i = 0; i < bandWidth.size(); i++)
    {
        std::cout<<"sd = "<<standardDeviation[i]<<"  "
                 <<"ssd = "<<sampleStandardDeviation[i]<<"  "
                 <<"bw = "<<bandWidth[i]<<"  "
                 <<"ssd/bw = "<<sampleStandardDeviation[i]/bandWidth[i]<<"  "
                 <<"silverman = "<<1.06*sampleStandardDeviation[i]/std::pow(numSamples, 1.0/5.0)<<"  "
                 <<"ratio = "<<1.06*sampleStandardDeviation[i]/std::pow(numSamples, 1.0/5.0)/bandWidth[i]<<"  "
                 <<std::endl;
    }*/
}

