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
#include <random>

#include <gtest/gtest.h>

#include "statistics/histogram_density_estimator.hpp"


/*!
 * \brief Test fixture for the HistogramDensityEstimator.
 */
class HistogramDensityEstimatorTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        HistogramDensityEstimatorTest()
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
 * This test checks that breakpoints are correctly constructed to cover the 
 * entire data range, be correctly ordered, and equidistant with spacing
 * equal to bin width. The test is carried out for a wide range of bin widths.
 */
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorBreaksTest)
{
    // floating point comparison tolerance:
    real eps = 10*std::numeric_limits<real>::epsilon();

    // need to sort test data manually in this test:
    std::sort(testData_.begin(), testData_.end());

    // create histogram estimator and set bin width:
    HistogramDensityEstimator hde;

    // some bin widths to try:
    std::vector<real> binWidth = {1.0, 0.1, 1e-5, std::sqrt(2.0)};

    // loop over bin widths and check that breaks are correct:
    for(auto it = binWidth.begin(); it != binWidth.end(); it++)
    {
        // set new bin width:
        hde.setBinWidth(*it);

        // create break points:
        std::vector<real> breaks = hde.createBreaks(
                testData_);

        // check that break points cover data range:
        ASSERT_LT(breaks.front(), 
                  *std::min_element(testData_.begin(), testData_.end()));
        ASSERT_GT(breaks.back(), 
                  *std::max_element(testData_.begin(), testData_.end()));

        // check that break points are increasing and  spacing is bin width:
        for(size_t i = 0; i < breaks.size() - 1; i++)
        {
            ASSERT_LT(breaks.at(i), breaks.at(i+1));
            ASSERT_NEAR(*it, breaks.at(i+1) - breaks.at(i), eps);
        }
    }
}


/*!
 * This test checks that the raw density estimate prior to spline interpolation
 * yields valid values. In particular, it ascertains that the density in each
 * bin is greater or equal to zero, that the sum of densities over all bins is
 * equal to one, and that the density in the first and last bin is zero. 
 * Additionally, it checks that exceptions are thrown if no valid bandwidth is
 * set.
 */
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorDensityTest)
{
    // floating point comparison tolerance:
    real eps = 10*std::numeric_limits<real>::epsilon();

    // need to sort test data manually in this test:
    std::sort(testData_.begin(), testData_.end());

    // create histogram estimator and set bin width:
    HistogramDensityEstimator hde;

    // check that exception is thrown if no bin width is set:
    ASSERT_THROW(hde.estimate(testData_), std::logic_error);

    // check that setting a negative of zero bin width does not work:
    ASSERT_THROW(hde.setBinWidth(-1.0), std::logic_error);
    ASSERT_THROW(hde.setBinWidth(0.0), std::logic_error);

    // some bin widths to try:
    // NOTE: if bin width too small, spline interpolation will fail!
    std::vector<real> binWidth = {1.0, 0.1, 1e-2, 0.1*std::sqrt(2.0)};

    // perform tests for each bin width value:
    for(auto bw : binWidth)
    {
        // set bin width:
        hde.setBinWidth(bw);

        // create breaks and estimate density:
        std::vector<real> breaks = hde.createBreaks(testData_);       
        std::vector<real> density = hde.calculateDensity(
                testData_,
                breaks);

        // ensure that density values are semi-positive definite:
        for(auto d : density)
        {
            ASSERT_LE(0.0, d);
        }

        // densities should sum to one:
        real sum = 0.0;
        for(auto d : density)
        {
            sum += d;
        }
        ASSERT_NEAR(1.0, sum, eps);

        // first and last bin should be empty by construction:
        ASSERT_NEAR(0.0, density.front(), 0.0);
        ASSERT_NEAR(0.0, density.back(), 0.0);
    }
}


/*!
 * Tests the behaviour in case of an empty input dataset, in which case an
 * all-zero density is returned.
 */
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorEmptyDatasetTest)
{
    // create empty dataset:
    std::vector<real> emptyData;

    // create histogram estimator and set bindwidth:
    HistogramDensityEstimator hde;

    // some bin widths to try:
    // NOTE: if bin width too small, spline interpolation will fail!
    std::vector<real> binWidth = {1.0, 0.1, 1e-2, 0.1*std::sqrt(2.0)};

    // perform tests for each binwidth value:
    for(auto bw : binWidth)
    {
        // set bin width:
        hde.setBinWidth(bw);

        // create breaks and estimate density:
        std::vector<real> breaks = hde.createBreaks(emptyData);       
        std::vector<real> density = hde.calculateDensity(
                emptyData,
                breaks);

        // ensure that density values are all zero:
        for(auto d : density)
        {
            ASSERT_EQ(0.0, d);
        }

        // first and last bin should be empty by constructrion:
        ASSERT_NEAR(0.0, density.front(), 0.0);
        ASSERT_NEAR(0.0, density.back(), 0.0);
    }
}


/*!
 * This test checks that the public interface for histogram density estimation 
 * yields a valid probability density. In particular, it makes sure that the
 * density outside the data range is zero, that the density within the data 
 * range is positive semi-definite, and that the integral over the density
 * spline equals one.
 */
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorEstimateTest)
{
    // floating point comparison tolerance:
    real eps = std::numeric_limits<real>::epsilon();

    // create histogram estimator and set bin width:
    HistogramDensityEstimator hde;

    // check that exception is thrown if no bin width is set:
    ASSERT_THROW(hde.estimate(testData_), std::logic_error);

    // check that setting a negative of zero bin width does not work:
    ASSERT_THROW(hde.setBinWidth(-1.0), std::logic_error);
    ASSERT_THROW(hde.setBinWidth(0.0), std::logic_error);

    // some bin widths to try:
    // NOTE: if bin width too small, spline interpolation will fail!
    std::vector<real> binWidth = {1.0, 0.1, 1e-2, 0.1*std::sqrt(2.0)};
    binWidth = {1.0, 0.1, 0.01, 0.1*std::sqrt(2.0)};

    // get data range:
    real dataMin = *std::min_element(testData_.begin(), testData_.end());
    real dataMax = *std::max_element(testData_.begin(), testData_.end());

    // check density computed for all values of bin width:
    for(auto bw : binWidth)
    {
        // estimate density at this bin width:
        hde.setBinWidth(bw);
        SplineCurve1D densitySpline = hde.estimate(testData_);

        // number of evaluation points:
        // (number is large so that probability integral quadrature is simple)
        int numEval = 10000;

        // evaluate density above data range:
        for(int i = 0; i < numEval; i++)
        {
            // eval point above data range:
            real eval = i*bw + dataMax + 2.0*bw;
            
            // evaluate density:
            real density = densitySpline.evaluate(
                    eval,
                    0);

            // this should always be zero:
            ASSERT_NEAR(0.0, density, std::sqrt(eps));
        }
 
        // evaluate density below data range:
        for(int i = 0; i < numEval; i++)
        {
            // eval point below data range: 
            real eval = -i*bw + dataMin - 0.5*bw;
            
            // evaluate density:
            real density = densitySpline.evaluate(
                    eval,
                    0);

            // this should always be zero:
            ASSERT_NEAR(0.0, density, std::sqrt(eps));
        }

        // evaluate spline within data range:
        std::vector<real> evalPoints;
        std::vector<real> evalDensities;
        real dataRange = dataMax - dataMin;
        real evalStep = (dataRange + 2.0*bw ) /(numEval);
        for(int i = 0; i < numEval; i++)
        {
            // calculate evaluation point:
            real eval = dataMin - bw + i*evalStep;
            evalPoints.push_back(eval);

            // evaluate density at this point:
            real density = densitySpline.evaluate(
                    eval,
                    0);
            evalDensities.push_back(density);
        }

        // assert positive semi-definiteness:
        for(auto d : evalDensities)
        {
            ASSERT_LE(0.0, d);
        }

        // check that integral over data range equals one:
        real integral = 0.0;
        for(auto d : evalDensities)
        {
            integral += d;
        }
        integral *= evalStep;
        ASSERT_NEAR(1.0, integral, std::sqrt(eps));
   }
}

