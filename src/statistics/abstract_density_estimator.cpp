#include "statistics/abstract_density_estimator.hpp"


/*!
 * Constructor sets all parameters to meaningless values and all flags to 
 * false.
 */
DensityEstimationParameters::DensityEstimationParameters()
    : binWidth_(-1.0)
    , binWidthIsSet_(false)
{

}


/*!
 * Sets the value of the bin width parameter to the given value and the 
 * corresponding flag to true.
 */
void
DensityEstimationParameters::setBinWidth(
        real binWidth)
{
    binWidth_ = binWidth;
    binWidthIsSet_ = true;
}


/*!
 * Returns the bin width parameter.
 */
real
DensityEstimationParameters::binWidth() const
{
    return binWidth_;
}


/*!
 * Returns a flag indicating whether the bin width parameter has been set.
 */
bool
DensityEstimationParameters::binWidthIsSet() const
{
    return binWidthIsSet_;
}

