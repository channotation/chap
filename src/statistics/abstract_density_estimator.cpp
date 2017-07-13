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
 * Sets the value of the band width parameter to the given value and the 
 * corresponding flag to true.
 */
void
DensityEstimationParameters::setBandWidth(
        real bandWidth)
{
    bandWidth_ = bandWidth;
    bandWidthIsSet_ = true;
}


/*!
 * Sets the value of the maximum permitted distance between two evaluation 
 * points.
 */
void
DensityEstimationParameters::setMaxEvalPointDist(
        real maxEvalPointDist)
{
    maxEvalPointDist_ = maxEvalPointDist;
}


/*!
 * Sets evaluation range cutoff to the given value. Note that this will by 
 * multiplied by the bandwidth.
 */
void
DensityEstimationParameters::setEvaluationRangeCutoff(
        real evalRangeCutoff)
{
    evalRangeCutoff_ = evalRangeCutoff;    
}


/*!
 * Sets the kernel function to the given value and the corresponding flag to 
 * true.
 */
void
DensityEstimationParameters::setKernelFunction(
        eKernelFunction kernelFunction)
{
    kernelFunction_ = kernelFunction;
    kernelFunctionIsSet_ = true;
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


/*!
 * Returns the band width parameter.
 */
real
DensityEstimationParameters::bandWidth() const
{
    return bandWidth_;
}


/*!
 * Returns a flag indicating whether the band width parameter has been set.
 */
bool
DensityEstimationParameters::bandWidthIsSet() const
{
    return bandWidthIsSet_;
}


/*!
 * Returns the maximum distance between two evaluation points.
 */
real
DensityEstimationParameters::maxEvalPointDist() const
{
    return maxEvalPointDist_;
}


/*!
 * Returns a flag indicating whether the maximum distance between two 
 * evaluation points has been set.
 */
bool
DensityEstimationParameters::maxEvalPointDistIsSet() const
{

}


/*!
 * Returns the value of the evaluation range cutoff.
 */
real
DensityEstimationParameters::evalRangeCutoff() const
{
    return evalRangeCutoff_;
}


/*!
 * Returns a flag indicating whether the evaluation range cutoff has been set.
 */
bool
DensityEstimationParameters::evalRangeCutoffIsSet() const
{
    return evluationRangeCutoffIsSet_;
}


/*!
 * Returns a flag indicating wether kernel function has been set.
 */
eKernelFunction
DensityEstimationParameters::kernelFunction() const
{
    return kernelFunction_;
}


/*!
 * Returns a function indicating whether kernel function has been set.
 */
bool
DensityEstimationParameters::kernelFunctionIsSet() const
{
    return kernelFunctionIsSet_;
}

