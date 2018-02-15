// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#include <random>

#include <gtest/gtest.h>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"

/*!
 * \brief Test fixture for the AmiseOptimalBandWidthEstimator.
 */
class AmiseOptimalBandWidthEstimatorTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        AmiseOptimalBandWidthEstimatorTest()
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
TEST_F(AmiseOptimalBandWidthEstimatorTest, 
       AmiseOptimalBandWidthEstimatorGaussianDensitySilvermanTest)
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
                    AmiseOptimalBandWidthEstimator bwe;
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

