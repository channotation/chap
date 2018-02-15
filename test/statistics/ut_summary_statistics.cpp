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
#include <cmath>
#include <numeric>
#include <limits>

#include <gtest/gtest.h>

#include "statistics/summary_statistics.hpp"


/*!
 * \brief Test fixture for SummaryStatistics.
 *
 * Initialises a simple hard-coded data set used in all tests.
 */
class SummaryStatisticsTest : public ::testing::Test
{
    public:

        // constructor for creating test data:
        SummaryStatisticsTest()
        {
            testData_ = {0.3, 1.5, -0.9, std::sqrt(2.0), -5.1};
        }

    
    protected:

        // test data:
        std::vector<real> testData_;
};


/*!
 * Checks that the minimum and maximum computed using SummaryStatistics agree
 * with the values obtained by directly searching the entire dataset.
 */
TEST_F(SummaryStatisticsTest, SummaryStatisticsMinMaxTest)
{
    // tolerance threshold for floating point comparison:
    real eps = std::numeric_limits<real>::epsilon();

    // create summary statistic object and update it with all available data:
    SummaryStatistics testDataSummary;
    for(size_t i = 0; i < testData_.size(); i++)
    {
        testDataSummary.update(testData_.at(i));
    }

    // manually compute minimum and maximum:
    real min = *(std::min_element(testData_.begin(), testData_.end()));
    real max = *(std::max_element(testData_.begin(), testData_.end()));

    // assert correctness:
    ASSERT_NEAR(min, testDataSummary.min(), eps);
    ASSERT_NEAR(max, testDataSummary.max(), eps);
}


/*!
 * Checks if the mean computed using SummaryStatistics agrees with the mean
 * computed by directly averaging the entire dataset.
 */
TEST_F(SummaryStatisticsTest, SummaryStatisticsMeanTest)
{
    // tolerance threshold for floating point comparison:
    real eps = std::numeric_limits<real>::epsilon();

    // create summary statistic object and update it with all available data:
    SummaryStatistics testDataSummary;
    for(size_t i = 0; i < testData_.size(); i++)
    {
        testDataSummary.update(testData_.at(i));
    }

    // manually compute mean:
    real mean = std::accumulate(testData_.begin(), testData_.end(), 0.0);
    mean /= testData_.size();

    // assert correctness:
    ASSERT_NEAR(mean, testDataSummary.mean(), eps);
}


/*!
 * Checks if the variance and standard deviation computed using 
 * SummaryStatistic agree with values obtained from calculating variance and
 * standard deviation directly from the overall data set.
 */
TEST_F(SummaryStatisticsTest, SummaryStatisticsVarSdTest)
{
    // tolerance threshold for floating point comparison:
    real eps = std::numeric_limits<real>::epsilon();

    // create summary statistic object and update it with all available data:
    SummaryStatistics testDataSummary;
    for(size_t i = 0; i < testData_.size(); i++)
    {
        testDataSummary.update(testData_.at(i));
    }

    // manually compute mean (needed to compute variance):
    real mean = std::accumulate(testData_.begin(), testData_.end(), 0.0);
    mean /= testData_.size();

    // manually compute variance and standard deviation:
    std::vector<real> diff(testData_.size());
    std::transform(
            testData_.begin(), 
            testData_.end(),
            diff.begin(),
            [mean](real x){ return x - mean; });
    real var = std::inner_product(
            diff.begin(), 
            diff.end(), 
            diff.begin(), 
            0.0);
    var /= testData_.size() - 1;
    real sd = std::sqrt(var);

    // assert correctness:
    ASSERT_NEAR(var, testDataSummary.var(), eps);
    ASSERT_NEAR(sd, testDataSummary.sd(), eps);
}

