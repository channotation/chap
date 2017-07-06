#include <algorithm>
#include <limits>

#include <gtest/gtest.h>

#include "statistics/summary_statistics.hpp"


/*!
 * Test fixture for SummaryStatistics. Initialises a simple hard-coded data set
 * used in all tests.
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

    // manually compute variance and standrd deviation:
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

