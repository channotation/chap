#include <algorithm>
#include <random>

#include <gtest/gtest.h>

#include "statistics/histogram_density_estimator.hpp"


/*
 *
 */
class HistogramDensityEstimatorTest : public ::testing::Test
{
    public:

        // constructor to set up test data:
        HistogramDensityEstimatorTest()
        {
            // parameters of normal distribution:
            real mu = -1.0;
            real sd = 3.0;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distribution(mu, sd);

            // create a random sample:
            int numSamples = 1e5;
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
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorBreaksTest)
{
    // floating point comparison tolerance:
    real eps = 10*std::numeric_limits<real>::epsilon();

    // need to sort test data manually in this test:
    std::sort(testData_.begin(), testData_.end());

    // create histogram estimator and set bindwidth:
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
                testData_.front(),
                testData_.back());

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


/*
 *
 */
TEST_F(HistogramDensityEstimatorTest, HistogramDensityEstimatorEstimateTest)
{
    // floating point comparison tolerance:
    real eps = 10*std::numeric_limits<real>::epsilon();


    // create histogram estimator and set bindwidth:
    HistogramDensityEstimator hde;

    // check that exception is thrown if no bin width is set:
    ASSERT_THROW(hde.estimate(testData_), std::logic_error);

    // check that setting a negative of zero bin width does not work:
    ASSERT_THROW(hde.setBinWidth(-1.0), std::logic_error);
    ASSERT_THROW(hde.setBinWidth(0.0), std::logic_error);

    // some bin widths to try:
    std::vector<real> binWidth = {1.0, 0.1, 1e-5, std::sqrt(2.0)};

    // loop over bin widths and check that breaks are correct:
    for(auto it = binWidth.begin(); it != binWidth.end(); it++)
    {
        // set bin width:
        hde.setBinWidth(*it);

        // estimate density:
        SplineCurve1D density = hde.estimate(testData_);

    }
}

