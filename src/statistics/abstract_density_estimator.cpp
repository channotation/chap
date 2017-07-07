#include "statistics/abstract_density_estimator.hpp"




/*
 *
 */
void
DensityEstimatorParameters::setBinWidth(
        real binWidth)
{
    binWidth_ = binWidth;
}


/*
 *
 */
real
DensityEstimatorParameters::binWidth() const
{
    return binWidth_;
}

