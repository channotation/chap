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


#include <algorithm>
#include <limits>
#include <random>

#include <gtest/gtest.h>

#include "statistics/gaussian_density_derivative.hpp"


/*!
 * \brief Test fixture for the GaussianDensityDerivative.
 */
class GaussianDensityDerivativeTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        GaussianDensityDerivativeTest()
        {
            // parameters of normal distribution:
            real muA = -1.0;
            real sdA = 0.5;
            real muB = 0.0;
            real sdB = 0.1;
            real muC = 5.0;
            real sdC = 1.5;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distributionA(muA, sdA);
            std::normal_distribution<real> distributionB(muB, sdB);
            std::normal_distribution<real> distributionC(muC, sdC);

            // create a random sample:
            size_t numSamples = 1e2 / 3;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distributionA(generator) );
                testData_.push_back( distributionB(generator) );
                testData_.push_back( distributionC(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
};


/*!
 * Checks that the shiftAndScale() functions maps correctly to the unit
 * interval and that the shiftAndScaleInverse() function restores the original
 * data points.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeShiftScaleTest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // set up data vectors:
    std::vector<real> vecA = {-1.0, 0.3, -0.215, 0.5, 1.0, 2.0};
    std::vector<real> vecB = {0.333, 0.891, 1.5, 10.0, 1.1, 2.7};

    // make copies for later reference:
    std::vector<real> refVecA = vecA;
    std::vector<real> refVecB = vecB;

    // find parameters for shifting:
    GaussianDensityDerivative gdd;
    auto ss = gdd.getShiftAndScaleParams(vecA, vecB);

    // map to unit interval:
    gdd.shiftAndScale(vecA, ss.first, ss.second);
    gdd.shiftAndScale(vecB, ss.first, ss.second);

    // check that all values are in unit interval:
    for(auto a : vecA)
    {
        ASSERT_GE(1.0, a);
        ASSERT_LE(0.0, a);
    }
    for(auto b : vecB)
    {
        ASSERT_GE(1.0, b);
        ASSERT_LE(0.0, b);        
    }

    // map back to original interval:
    gdd.shiftAndScaleInverse(vecA, ss.first, ss.second);
    gdd.shiftAndScaleInverse(vecB, ss.first, ss.second);
    
    // check that values match original values:
    for(size_t i = 0; i < vecA.size(); i++)
    {
        ASSERT_NEAR(vecA[i], refVecA[i], tolerance);
    }
    for(size_t i = 0; i < vecB.size(); i++)
    {
        ASSERT_NEAR(vecB[i], refVecB[i], tolerance);
    }
}


/*!
 * Checks that space partitioning produces correct centres and that each data
 * point is associated with the correct centre.
 */
TEST_F(GaussianDensityDerivativeTest, 
       GaussianDensityDerivativeSpacePartitioningTest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // run this test for a variaty of bandwidths:
    std::vector<real> bandwidth = {2.0, 1.0, 0.5, 0.3, 0.1, 0.0001};
    for(auto bw : bandwidth)
    {
        // sample data in interval [0,1]:
        std::vector<real> sample = {0.0, 0.33, 0.5, 0.7, 0.4, 0.5, 0.121, 0.9, 1.0};
    
        // prepare density derivative estimator:
        GaussianDensityDerivative gdd;
        gdd.setDerivOrder(2);
        gdd.setBandWidth(bw);
        gdd.setErrorBound(0.001);

        gdd.q_ = gdd.setupCoefQ(sample.size());
        gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
        gdd.rc_ = gdd.setupCutoffRadius();
        gdd.trunc_ = gdd.setupTruncationNumber();

        // obtain cluster centres:
        std::vector<real> centres = gdd.setupClusterCentres();
        gdd.centres_ = centres;

        // largest and smallest cluster centre:
        real cMin = *std::min_element(centres.begin(), centres.end());
        real cMax = *std::max_element(centres.begin(), centres.end());

        // check that all centres fall into unit interval:
        ASSERT_GE(cMax, cMin);
        ASSERT_GE(cMin, 0.0);
        ASSERT_GE(1.0, cMax);

        // check that centre distance is bounded by half bandwidth: 
        for(size_t i = 0; i < centres.size() - 1; i++)
        {
            ASSERT_LE(centres[i+1] - centres[i], bw/2.0 + tolerance);
        }

        // find interval centre for each data point:
        std::vector<unsigned int> idx = gdd.setupClusterIndices(sample);

        // for each data point check that correct centre was found: 
        for(size_t i = 0; i < sample.size(); i++)
        {
            // find distance from closest centre on record:
            real d = sample.at(i) - centres.at(idx.at(i));
            
            // calculate distance from sample to all centres:
            std::vector<real> dist;
            dist.reserve(centres.size());
            for(auto c : centres)
            {
                dist.push_back(std::abs(sample[i] - c));
            }

            // find minimal distance:
            auto minIter = min_element(dist.begin(), dist.end());
            real minDist = *minIter;

            // check that centre on record is actually closest:
            ASSERT_LE(d, minDist);
        }
    }
}


/*!
 * Checks the the correct truncation number is found for a variety of 
 * parameters and random data sets.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeTruncationTest)
{
    // parameters of normal distribution:
    real muA = -1.0;
    real sdA = 0.5;
    real muB = 0.0;
    real sdB = 0.1;
    real muC = 5.0;
    real sdC = 1.5;

    // prepare random distribution:
    std::default_random_engine generator;
    std::normal_distribution<real> distributionA(muA, sdA);
    std::normal_distribution<real> distributionB(muB, sdB);
    std::normal_distribution<real> distributionC(muC, sdC);
    
    // perform test for various parameters:
    std::vector<unsigned int> numSamples = {10, 100, 1000, 10000};
    std::vector<real> bandwidth = {1.0, 0.1, 0.01, 0.001, 0.0001};
    std::vector<real> epsilon = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
    for(auto n : numSamples)
    {
        // create a random sample:
        std::vector<real> sample;
        for(size_t i = 0; i < n; i++)
        {
            sample.push_back( distributionA(generator) );
            sample.push_back( distributionB(generator) );
            sample.push_back( distributionC(generator) );
        }
        
        // loop over parameters:
        for(auto bw : bandwidth)
        {
            for(auto eps : epsilon)
            {
                // set parameters:
                GaussianDensityDerivative gdd;
                unsigned int deriv = 2;
                gdd.setDerivOrder(deriv);
                gdd.setBandWidth(bw);
                gdd.setErrorBound(eps);

                // determine required coefficients:
                gdd.q_ = gdd.setupCoefQ(sample.size());
                gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
                gdd.rc_ = gdd.setupCutoffRadius();
                gdd.centres_ = gdd.setupClusterCentres();

                // calculate truncation number:
                unsigned int trunc = gdd.setupTruncationNumber();

                // check validity of truncation number:
                ASSERT_FALSE(std::isnan(trunc));
                ASSERT_GE(trunc, 0);

                // check that truncation criterion is met:
                real b = std::min<real>(gdd.rc_, 
                                  0.5*(gdd.ri_ + std::sqrt(gdd.ri_*gdd.ri_ + 
                                  8.0*trunc*gdd.bw_*gdd.bw_)));
                real error = std::sqrt(gdd.factorial(gdd.r_))
                           / gdd.factorial(trunc)
                           * std::pow((gdd.ri_ * b / (gdd.bw_*gdd.bw_)), trunc)
                           * std::exp(-(std::pow(gdd.ri_ - b, 2))
                                      / (gdd.bw_*gdd.bw_));
                ASSERT_LE(error, gdd.epsPrime_);                
            }
        }
    }
}


/*!
 * Checks that the a-coefficients are correctly computed for two different
 * derivative orders by comparison to manually calculated coefficients.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeCoefATest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // prepare density derivative estimator:
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(0.001);

    // manually computed coefficients for this parameter vector:
    std::vector<real> coefATrue = {1, -2, 1, -1};

    // calculate coefficients using reference and optimised implementation:
    std::vector<real> coefA = gdd.setupCoefA();

    // check correct number of coefficients:
    ASSERT_EQ(coefATrue.size(), coefA.size());

    // check correct value of coefficients:
    for(size_t i = 0; i < coefATrue.size(); i++)
    {
        ASSERT_NEAR(coefATrue[i], coefA[i], tolerance);
    }

    // set a different derivative order:
    gdd.setDerivOrder(5);

    // manually computed coefficients for this parameter vector:
    coefATrue = {  1, -5,  10, -10,  5,  -1, 
                 -10, 30, -30,  10, 15, -15};

    // calculate coefficients using reference and optimised implementation:
    coefA = gdd.setupCoefA();
    
    // check correct number of coefficients:
    ASSERT_EQ(coefATrue.size(), coefA.size());

    // check correct value of coefficients:
    for(size_t i = 0; i < coefATrue.size(); i++)
    {
        ASSERT_NEAR(coefATrue[i], coefA[i], tolerance);
    }
}


/*!
 * Checks that all B coefficients computed for a variety of parameters and 
 * random data sets of various sizes are finate and not NaN or Inf. This 
 * might happen for large truncation numbers where floating point numbers might
 * over or underflow.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeCoefBTest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // parameters of normal distribution:
    real muA = -1.0;
    real sdA = 0.5;
    real muB = 0.0;
    real sdB = 0.1;
    real muC = 5.0;
    real sdC = 1.5;

    // prepare random distribution:
    std::default_random_engine generator;
    std::normal_distribution<real> distributionA(muA, sdA);
    std::normal_distribution<real> distributionB(muB, sdB);
    std::normal_distribution<real> distributionC(muC, sdC);
 
    // perform test for various parameters:
    std::vector<unsigned int> numSamples = {10, 100, 1000, 10000};
    std::vector<real> bandwidth = {1.0, 0.1, 0.01, 0.001, 0.0001};
    std::vector<real> epsilon = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
    for(auto n : numSamples)
    {
        // create a random sample:
        std::vector<real> sample;
        for(size_t i = 0; i < n; i++)
        {
            sample.push_back( distributionA(generator) );
            sample.push_back( distributionB(generator) );
            sample.push_back( distributionC(generator) );
        }
       
        // shift and scale data:
        GaussianDensityDerivative gdd;
        std::vector<real> eval = sample; // as dummy input to scaling
        auto ss = gdd.getShiftAndScaleParams(sample, eval);
        gdd.shiftAndScale(sample, ss.first, ss.second);
 
        // loop over parameters:
        for(auto bw : bandwidth)
        {
            for(auto eps : epsilon)
            {
                // set parameters:
                unsigned int deriv = 2;
                gdd.setDerivOrder(deriv);
                gdd.setBandWidth(bw);
                gdd.setErrorBound(eps);

                // determine required coefficients:
                gdd.q_ = gdd.setupCoefQ(sample.size());
                gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
                gdd.rc_ = gdd.setupCutoffRadius();
                gdd.centres_ = gdd.setupClusterCentres();
                gdd.idx_ = gdd.setupClusterIndices(sample);
                gdd.trunc_ = gdd.setupTruncationNumber();

                // compute coefficients:
                std::vector<real> coefB = gdd.setupCoefB(sample);
                
                // check validity of coefficients:
                for(unsigned int i = 0; i < coefB.size(); i++)
                {
                    ASSERT_FALSE(std::isnan(coefB[i]));
                    ASSERT_FALSE(std::isinf(coefB[i]));
                }
            }
        }
    }
}


/*!
 * Checks that the derivative computed using the approximate method of Raykar 
 * is consistent with the derivative computed by direct evaluation of the 
 * double sum. This is done for various parameters combinations and for a 
 * number of different size random data sets sampled from a Gaussian mixture
 * distribution.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeConsistencyTest)
{
    // parameters of normal distribution:
    real muA = -1.0;
    real sdA = 0.5;
    real muB = 0.0;
    real sdB = 0.1;
    real muC = 5.0;
    real sdC = 1.5;

    // prepare random distribution:
    std::default_random_engine generator;
    std::normal_distribution<real> distributionA(muA, sdA);
    std::normal_distribution<real> distributionB(muB, sdB);
    std::normal_distribution<real> distributionC(muC, sdC);
    
    // carry out test on various sized random samples:
    // can not be much larger as direct evaluation too slow:
    std::vector<real> numSamples = {5, 50, 500};
    for(auto n : numSamples)
    {
        // create a random sample:
        std::vector<real> sample;
        for(size_t i = 0; i < n; i++)
        {
            sample.push_back( distributionA(generator) );
            sample.push_back( distributionB(generator) );
            sample.push_back( distributionC(generator) );
        }

        // evaluation points equal to sample points:
        std::vector<real> eval = sample;
        
        // map input data to unit interval:
        GaussianDensityDerivative gdd;
        auto ss = gdd.getShiftAndScaleParams(sample, eval);
        gdd.shiftAndScale(sample, ss.first, ss.second);
        gdd.shiftAndScale(eval, ss.first, ss.second);

        // carry out test for multiple parameter combinations:
        std::vector<real> epsilon = {1e-1, 1e-2, 1e-3, 1e-4};
        std::vector<real> bandwidth = {10.0, 1.0, 0.1, 0.01, 0.001};
        for(auto eps : epsilon)
        {
            for(auto bw : bandwidth)
            {
                // set parameters:
                gdd.setDerivOrder(2);
                gdd.setBandWidth(bw);
                gdd.setErrorBound(eps);

                // estimate derivative via direct loop and via approximate method:
                std::vector<real> derivDirect = gdd.estimateDirect(sample, eval);
                std::vector<real> derivApprox = gdd.estimateApprox(sample, eval);

                // check equality of estimation methods:
                for(size_t i = 0; i < eval.size(); i++)
                {
                    // difference between estimation methods:
                    real tol = std::max(eps,
                                        eps*std::fabs(derivDirect[i]));

                    // assertion on relative error:
                    ASSERT_NEAR(derivDirect[i], derivApprox[i], 5*tol);
                }
            }
        }
    }    
}

