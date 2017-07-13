#include <algorithm>
#include <cmath>

#include "geometry/linear_spline_interp_1D.hpp"
#include "statistics/kernel_density_estimator.hpp"


/*
 *
 */
SplineCurve1D
KernelDensityEstimator::estimate(
        const std::vector<real> &samples)
{
    // construct evaluation points:
    std::vector<real> evalPoints = createEvaluationPoints(
            samples);
   
    // determine density at evaluation points:
    std::vector<real> density = calculateDensity(
            samples,
            evalPoints);

    // interpolate density between sample points:
    LinearSplineInterp1D Interp;
    return Interp(evalPoints, density);
}


/*
 *
 */
void
KernelDensityEstimator::setParameters(
        const DensityEstimationParameters &params)
{
    // check if required parameters are set:
    if( params.kernelFunctionIsSet() )
    {
        setKernelFunction(params.kernelFunction());
    }
    else
    {
        throw std::runtime_error("Kernel function is not set!");
    } 

    if( params.bandWidthIsSet() )
    {
        setBandWidth(params.bandWidth());
    }
    else
    {
        throw std::runtime_error("Kernel bandwidth is not set!");
    }

    if( params.evalRangeCutoffIsSet() )
    {
        setEvalRangeCutoff(params.evalRangeCutoff());
    }
    else
    {
        throw std::runtime_error("Evaluation range cutoff is not set!");
    }

    if( params.maxEvalPointDistIsSet() )
    {
        setMaxEvalPointDist(params.maxEvalPointDist());
    }
    else
    {
        throw std::runtime_error("Maximum evluation point distance is not set!");
    }
}


/*!
 * Sets the band width to the given value. Throws exception if band width is
 * not positive.
 */
void 
KernelDensityEstimator::setBandWidth(
        const real bandWidth)
{
    // sanity check
    if( bandWidth <= 0.0 )
    {
        throw std::logic_error("Bandwidth may not be negative.");
    }

    // set internal band width parameter:
    bandWidth_ = bandWidth;
}


/*!
 * Sets the maximum evaluation point distance parameter to the given value. 
 * Throws an exception if value is not positive.
 */
void
KernelDensityEstimator::setMaxEvalPointDist(
        const real maxEvalPointDist)
{
    // sanity check:
    if( maxEvalPointDist <= 0.0 )
    {
        throw std::logic_error("Maximum distance between evaluation points "
        "must be positive.");
    }

    // set avlue of internal parameter:
    maxEvalPointDist_ = maxEvalPointDist;
}


/*!
 * Sets the evaluation range cutoff to a given value. Throws an exception if
 * the value is negative.
 */
void
KernelDensityEstimator::setEvalRangeCutoff(
        const real evalRangeCutoff)
{
    // sanity check:
    if( evalRangeCutoff < 0 )
    {
        throw std::logic_error("Evaluation range cutoff may not be negative!");
    }

    // set internal parameter:
    evalRangeCutoff_ = evalRangeCutoff;
}


/*!
 * Sets kernel function to the given value.
 */
void
KernelDensityEstimator::setKernelFunction(
        const eKernelFunction kernelFunction)
{
    kernelFunction_ = kernelFunction;
}


/*!
 * Auxiliary function that creates a set of equidistant evaluation points at 
 * which the density will be evaluated. 
 *
 * The points are chosen such that the  entire data range plus an extra margin 
 * is covered. The extra margin is chosen as some factor times the bandwidth, 
 * where the factor can be set via setEvalRangeCutoff(). The margin is added
 * to each end of the data range. This is done so that the density at the 
 * endpoints decays to almost zero (for Kernels with local support, it can be 
 * exactly zero).
 *
 * The spacing of the evaluation points can be controlled using 
 * setMaxEvalPointDist(), which sets the upper limit for the distance between
 * two subsequent evaluation points. The real distance is calculated by 
 * requiring that the number of evaluation points be a power of two, in order
 * to facilitate the evaluation of the density using an FFT-based algorithm.
 */
std::vector<real>
KernelDensityEstimator::createEvaluationPoints(
        const std::vector<real> &samples)
{
    // find range covered by data:
    real rangeLo = *std::min_element(samples.begin(), samples.end());
    real rangeHi = *std::max_element(samples.begin(), samples.end());

    // extend this by multiple of bandwidth:
    rangeLo -= evalRangeCutoff_ * bandWidth_;
    rangeLo += evalRangeCutoff_ * bandWidth_;

    // find required number of evaluation points and corresponding step:
    real range = rangeHi - rangeLo;
    size_t numEvalPoints = calculateNumEvalPoints(range);
    real deltaEvalPoints = range/(numEvalPoints - 1);

    // create set of evalPoints:
    std::vector<real> evalPoints(numEvalPoints, rangeLo);
    for(size_t i = 0; i < numEvalPoints; i++)
    {
        evalPoints[i] += i*deltaEvalPoints;
    }

    // return evaluation points:
    return evalPoints;
}


/*!
 * Auxiliary function for calculating the number of evaluation points required
 * for covering a given range of \f$ L \f$ with equidistant evaluation points 
 * no further than a specified distance, \f$ \Delta L \f$, apart.
 *
 * The minimum number of points required to achieve this is calculated as
 * \f$ \lceil L/\Delta L \rceil + 1\f$. This number is is then rounded up to 
 * the nearest power of two. This is done to facilitate an FFT-based 
 * implementation of the kernel density estimator.
 */
size_t
KernelDensityEstimator::calculateNumEvalPoints(
        const real range)
{
    // get minimal number of points required to fully cover given range:
    size_t numPoints = std::ceil(range / maxEvalPointDist_) + 1;

    // round this up to the nearest power of two:
    std::pow(2, std::ceil(std::log2(numPoints)));

    // return number of points:
    return numPoints;
}


/*
 *
 */
std::vector<real>
KernelDensityEstimator::calculateDensity(
        const std::vector<real> &samples,
        const std::vector<real> &evalPoints)
{
    // create kernel:
    KernelFunctionPointer Kernel = KernelFunctionFactory::create(
            kernelFunction_);

    // initialise density vector as zero:
    std::vector<real> density(evalPoints.size(), 0.0);

    // normalisation constant:
    real normalisation = 1.0 / samples.size() * bandWidth_;
    normalisation *= Kernel -> normalisingFactor();

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints.size(); i++)
    {
        // density is sum over kernel distances:
        for(auto sample : samples)
        {
            density[i] += Kernel -> operator()( evalPoints[i] -sample );
        }

        // normalise density at this evaluation point:
        density[i] *= normalisation;
    }

    // return density:
    return(density);
}

