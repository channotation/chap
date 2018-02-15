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

#include "statistics/kernel_density_estimator.hpp"


/*!
 * \brief Test fixture for testing the KernelDensityEstimator. 
 *
 * Provides a simple vector of test data sampled from a Gaussian distribution.
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
            real mu = mu_;
            real sd = sd_;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distribution(mu, sd);

            // create a random sample:
            size_t numSamples = 100;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
        real sd_ = 0.1;
        real mu_ = -std::sqrt(2.0);
};


/*!
 * This test checks that exceptions are thrown in case of unset or nonsensical
 * parameters.
 */
TEST_F(KernelDensityEstimatorTest, KernelDensityEstimatorParameterTest)
{
    // create density estimator and parameter container:
    KernelDensityEstimator kde;
    DensityEstimationParameters params;

    // assert exception in case of unset parameters:
    ASSERT_THROW(kde.estimate(testData_), std::runtime_error);

    // assert exceptions parameters:
    ASSERT_THROW(kde.setMaxEvalPointDist(-1.0), std::logic_error);
    ASSERT_THROW(kde.setMaxEvalPointDist(0.0), std::logic_error);
    ASSERT_THROW(kde.setEvalRangeCutoff(-1.0), std::logic_error);

    // assert exceptions on unset parameters:
    ASSERT_THROW(kde.setParameters(params), std::runtime_error);
    params.setKernelFunction(eKernelFunctionGaussian);
    ASSERT_THROW(kde.setParameters(params), std::runtime_error);
    params.setBandWidth(1.0);
    ASSERT_THROW(kde.setParameters(params), std::runtime_error);
    params.setEvalRangeCutoff(3.0);
    ASSERT_THROW(kde.setParameters(params), std::runtime_error);
    params.setMaxEvalPointDist(1.0);
    ASSERT_NO_THROW(kde.setParameters(params));
}


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
        ASSERT_GE(
                dataLo + std::numeric_limits<real>::epsilon(), 
                evalPoints.front());
        ASSERT_LE(
                dataHi - std::numeric_limits<real>::epsilon(),
                evalPoints.back());
    }
}


/*!
 * This test checks that the kernel density estimator with a 
 * GaussianKernelFunction estimates a valid probability density. To be precise,
 * it asserts that all densities are non-negative. it does NOT check whether
 * the probability density integrates to one, as this is done by in the 
 * subsequent test on the interpolated density (this has the advantage that the
 * integration step can be chosen independent of the evaluation step and band
 * width).
 */
TEST_F(
        KernelDensityEstimatorTest, 
        KernelDensityEstimatorGaussianRawDensityTest)
{    
    // define a number of bandwidths to test:    
    std::vector<real> bandWidths = {10.0, 1.0, 0.1, 0.01, 0.001};

    // perform tests for all the above bandwidths:
    for(auto bw : bandWidths)
    {
        // create kernel density estimator and set parameters:
        // (note that large cutoff makes probability mass outside cutoff range
        // negligible).
        KernelDensityEstimator kde;
        kde.setBandWidth(bw);
        kde.setEvalRangeCutoff(3.0); 
        kde.setMaxEvalPointDist(sd_);
        kde.setKernelFunction(eKernelFunctionGaussian);

        // get evaluation points and calculate density:
        std::vector<real> evalPoints = kde.createEvaluationPoints(testData_);
        std::vector<real> density = kde.calculateDensity(testData_, evalPoints);

        // assert non-negativity of density:
        for(auto d : density)
        {
            ASSERT_LE(0.0, d);
        }
    }
}


/*!
 * Test the behaviour in the case of an empty input data set. In this case it 
 * is asserted that the estimated density is zero everywhere.
 */
TEST_F(
        KernelDensityEstimatorTest, 
        KernelDensityEstimatorEmptyDatasetTest)
{    
    // define a number of bandwidths to test:    
    std::vector<real> bandWidths = {10.0, 1.0, 0.1, 0.01, 0.001};

    // create empty data set:
    std::vector<real> emptyData;

    // perform tests for all the above bandwidths:
    for(auto bw : bandWidths)
    {
        // create kernel density estimator and set parameters:
        // (note that large cutoff makes probability mass outside cutoff range
        // negligible).
        KernelDensityEstimator kde;
        kde.setBandWidth(bw);
        kde.setEvalRangeCutoff(0.0); 
        kde.setMaxEvalPointDist(sd_);
        kde.setKernelFunction(eKernelFunctionGaussian);

        // get evaluation points and calculate density:
        std::vector<real> evalPoints = kde.createEvaluationPoints(emptyData);
        std::vector<real> density = kde.calculateDensity(emptyData, evalPoints);

        // assert non-negativity of density:
        for(auto d : density)
        {
            ASSERT_EQ(0.0, d);
        }
    }
}


/*! 
 * This test checks that the spline curve returned by the 
 * KernelDensityEstimator using a GaussianKernelFunction is a valid probability 
 * density. In particular, it asserts that it only evaluates to positive 
 * values, that it integrates to one and that it is zero outside of the data
 * range.
 */
TEST_F(
        KernelDensityEstimatorTest, 
        KernelDensityEstimatorGaussianInterpDensityTest)
{
    // set parameter ranges:
    std::vector<real> bandWidths = {10.0, 1.0, 1e-1, 1e-2, 1e-3};
    std::vector<real> evalPointDistanceFactors = {1.0, 1e-1, 1e-2};

    // conduct test for all bandwidths:
    for(auto bw : bandWidths)
    {
        // vary evaluation step relative to bandwidth:
        for(auto evalPointDistFac : evalPointDistanceFactors)
        {
            // create kernel density estimator and set parameters:
            KernelDensityEstimator kde;
            DensityEstimationParameters params;
            params.setBandWidth(bw);
            params.setBandWidthScale(1.0);
            params.setEvalRangeCutoff(5.0); 
            params.setMaxEvalPointDist(evalPointDistFac*bw);
            params.setKernelFunction(eKernelFunctionGaussian);
            kde.setParameters(params);

            // obtain density as spline curve:
            SplineCurve1D densitySpline = kde.estimate(testData_);

            // find density range:
            real margin = 3*bw;
            real rangeLo = densitySpline.knotVector().front() - margin;
            real rangeHi = densitySpline.knotVector().back() + margin;
            
            // create evaluation points covering data range:
            size_t numEvalPoints = 10000;
            real evalStep = (rangeHi - rangeLo) / (numEvalPoints - 1);
            std::vector<real> evalPoints;
            for(size_t i = 0; i < numEvalPoints; i++ )
            {
                evalPoints.push_back(rangeLo + i*evalStep);
            }
            
            // assert non-negativity of interpolated density:
            for(auto eval : evalPoints)
            {
                ASSERT_LE(0.0, densitySpline.evaluate(eval, 0));
            }
           
            // assert that density outside data range is zero:
            ASSERT_NEAR(
                    0.0, 
                    densitySpline.evaluate(rangeLo, 0),
                    std::numeric_limits<real>::epsilon());
            ASSERT_NEAR(
                    0.0, 
                    densitySpline.evaluate(rangeHi, 0),
                    std::numeric_limits<real>::epsilon());

            // assert that density integrates to one:
            real integral = 0.0;
            for(auto eval : evalPoints)
            {
                integral += densitySpline.evaluate(eval, 0);
            }
            integral *= evalStep;
            ASSERT_NEAR(
                    1.0, 
                    integral, 
                    std::sqrt(std::numeric_limits<real>::epsilon()));
        }
    }
}

