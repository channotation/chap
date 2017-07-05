#include <cmath>
#include <limits>

#include "statistics-utilities/summary_statistics.hpp"


/*!
 * Initialises all summary statistics and the sample counter as
 * zero.
 */
SummaryStatistics::SummaryStatistics()
    : min_(std::numeric_limits<real>::infinity())
    , max_(-std::numeric_limits<real>::infinity())
    , mean_(0.0)
    , sumSquaredMeanDiff_(0.0)
    , num_(0.0)
{
    
}


/*!
 * This method will update all summary statistics with the given new value and
 * also increment the sample counter.
 */
void
SummaryStatistics::update(const real newValue)
{
    // increment number of samples:
    num_++;

    // updating min and max is trivial:
    if( newValue > max_ )
    {
        max_ = newValue;
    }
    if( newValue < min_ )
    {
        min_ = newValue;
    }

    // update mean:
    real delta = newValue - mean_;
    mean_ += delta/num_;

    // update squared difference from mean:
    sumSquaredMeanDiff_ += delta*(newValue - mean_);
}


/*!
 * Getter method for obtaining the minimum value.
 */
real
SummaryStatistics::min() const
{
    return min_;
}


/*! 
 * Getter method for obtaining the maximum value.
 */
real
SummaryStatistics::max() const
{
    return max_;
}


/*!
 * Getter method for obtaining the (arithmetic) mean value.
 */
real
SummaryStatistics::mean() const
{
    return mean_;
}


/*!
 * Getter method for obtaining the variance. Note that this will return zero
 * if the estimate consists of less than two samples.
 */
real
SummaryStatistics::var() const
{
    if( num_ < 2 )
    {
        return 0.0;
    }
    else
    {
        return varFromSumSquaredMeanDiff();
    }
}


/*!
 * Getter method for obtaining the standard deviation.
 */
real
SummaryStatistics::sd() const
{
    if( num_ < 2 )
    {
        return 0.0;
    }
    else
    {
        return std::sqrt( varFromSumSquaredMeanDiff() );
    }
}


/*!
 * Getter method for obtaining the number of samples, i.e. the number of times
 * update() has been called.
 */
int
SummaryStatistics::num() const
{
    return num_;
}


/*!
 * Convenience function for converting sum of squared differences from mean to
 * variance. This is written as a separate function to be used with both the
 * standard deviation and variance getter methods and is inline for
 * performance.
 */
inline real
SummaryStatistics::varFromSumSquaredMeanDiff() const
{
    return sumSquaredMeanDiff_ / (num_ - 1.0);
}

