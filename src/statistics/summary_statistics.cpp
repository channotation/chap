#include <cmath>
#include <limits>
#include <stdexcept>

#include "statistics/summary_statistics.hpp"


/*!
 * Initialises all summary statistics and the sample counter as
 * zero. Minimum and maximum are initialised as the largest and smallest 
 * representable real number respectively. This is done instead of using 
 * infinities because JSON can not represent infinities.
 */
SummaryStatistics::SummaryStatistics()
    : min_(std::numeric_limits<real>::max())
    , max_(-std::numeric_limits<real>::max())
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
    // handle infinities:
    if( std::isinf(newValue) )
    {
        return;
    }

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
 * Convenience function to update a vector of SummaryStatistics with a vector
 * of new values.
 */
void
SummaryStatistics::updateMultiple(
        std::vector<SummaryStatistics> &stat,
        const std::vector<real> &newValues)
{
    // sanity check:
    if( stat.size() != newValues.size() )
    {
        throw std::logic_error("Can not update summary statistics vector with "
                               "data vector of different size.");
    }

    // update each value individually:
    for(size_t i = 0; i < stat.size(); i++)
    {
        stat[i].update(newValues[i]);
    }
}


/*!
 * Shifts the value of minimum, maximum, and mean by the given amount. Standard
 * deviation, variance, and number of samples are unaffacted. This is useful if
 * SummaryStatistics is used as a data container, but once shift() has been 
 * called, update() should no longer be called.
 */
void
SummaryStatistics::shift(
        const real shift)
{
    min_ += shift;
    max_ += shift;
    mean_ += shift;
}


/*!
 * Getter method for obtaining the minimum value.
 */
real
SummaryStatistics::min() const
{
    if( num_ > 0 )
    {
        return mendInfinity( min_ );
    }
    else
    {
        return -std::numeric_limits<real>::max();
    }
}


/*! 
 * Getter method for obtaining the maximum value.
 */
real
SummaryStatistics::max() const
{
    // handle case of no data:
    if( num_ > 0 )
    {
       return mendInfinity( max_ );
    }
    else
    {
        return std::numeric_limits<real>::max();
    }
}


/*!
 * Getter method for obtaining the (arithmetic) mean value.
 */
real
SummaryStatistics::mean() const
{
    return mendInfinity( mean_ );
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
        return mendInfinity( varFromSumSquaredMeanDiff() );
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
        return mendInfinity( std::sqrt( varFromSumSquaredMeanDiff() ) );
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


/*!
 * Auxiliary function for handling infinities. The getter methods for all
 * summary statistics call this function on the returned value and this 
 * functions turns negative and positive infinities into the smallest and 
 * largest real number respectively.
 */
real
SummaryStatistics::mendInfinity(real value) const
{
    if( std::isinf(value) )
    {
        if( value < 0.0 )
        {
            return -std::numeric_limits<real>::max();
        }
        else
        {
            return std::numeric_limits<real>::max();
        }
    }
    else
    {
        return value;
    }
}

