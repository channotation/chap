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

#include "geometry/linear_spline_interp_1D.hpp"
#include "statistics/kernel_density_estimator.hpp"


/*!
 * Public interface for density estimation. Returns a SplineCurve1D 
 * representing the kernel density estimate calculated for the given sample 
 * vector. Will throw an exception if called before parameters have been set
 * using setParameters().
 *
 * Internally, this will create a set of evaluation points covering the entire
 * data range plus some user specified margin using createEvaluationPoints().
 * It then uses a user specified AbstractKernelFunction to estimate the 
 * probability density of the samples and adds one additional density point at
 * either end of the evaluation range to ensure that the density is zero 
 * outside this range. The density points are then linearly interpolated and 
 * the result is returned as a SplineCurve1D.
 */
SplineCurve1D
KernelDensityEstimator::estimate(
        std::vector<real> &samples)
{
    // parameters set?
    if( !parametersSet_ )
    {
        throw std::runtime_error("kernel density estimation parameters not set!");
    }
 
    // construct evaluation points:
    std::vector<real> evalPoints = createEvaluationPoints(
            samples);
   
    // determine density at evaluation points:
    std::vector<real> density = calculateDensity(
            samples,
            evalPoints);

    // set endpoint density to zero:
    endpointDensityToZero(
            density,
            evalPoints);

    // interpolate density between sample points:
    LinearSplineInterp1D Interp;
    return Interp(evalPoints, density);
}


/*!
 * Implements the parameter setting interface defined in 
 * AbstractDensityEstimator. This function can be used to set the parameters
 * controlling the kernel density estimator. It will throw exceptions if any
 * required parameter is not set in the given DensityEstimationParameters
 * object or if any of the given parameters has a nonsensical value.
 *
 * In particular the following parameters are expected to be set:
 *
 * @param params.kernelFunction_ - an eKernelFunction specifying the kernel to 
 * be used
 * @param params.bandWidth_ - a real giving the kernel band width
 * @param params.evalRangeCutoff_ - a real specifying the cutoff of the 
 * evaluation range in multiples of the bandWidth
 * @param params.maxEvalPointDist - a real specifying the maximum distance 
 * between two subsequent evaluation points
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

    if( params.bandWidthScaleIsSet() )
    {
        setBandWidthScale(params.bandWidthScale());
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

    // set flag:
    parametersSet_ = true;
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
        throw std::logic_error("Bandwidth must be positive!");
    }

    // set internal band width parameter:
    bandWidth_ = bandWidth;
}


/*!
 * Sets the band width scake to the given value. Throws exception if scale is
 * not positive.
 */
void 
KernelDensityEstimator::setBandWidthScale(
        const real scale)
{
    // sanity check
    if( scale <= 0.0 )
    {
        throw std::logic_error("Band width scale must be positive!");
    }

    // set internal band width parameter:
    bandWidthScale_ = scale;
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
        "must be positive!");
    }

    // set value of internal parameter:
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
 *
 * This function also enforces a minimum of 512 sampling points to deal with 
 * situations where there are only very few data points very close to one 
 * another.
 */
std::vector<real>
KernelDensityEstimator::createEvaluationPoints(
        const std::vector<real> &samples)
{
    // handle special case of empty sample:
    real rangeLo = 0.0;
    real rangeHi = 0.0;
    if( samples.size() != 0 )
    {
        // find range covered by data:
        rangeLo = *std::min_element(samples.begin(), samples.end());
        rangeHi = *std::max_element(samples.begin(), samples.end());
    }

    // extend this by multiple of bandwidth:
    rangeLo -= evalRangeCutoff_ * bandWidth_;
    rangeHi += evalRangeCutoff_ * bandWidth_;

    // calculate data range:
    // (this enforces a minimum of 512 evaluation points)
    size_t minNumEvalPoints = 512;
    real range = rangeHi - rangeLo;
    if( range < minNumEvalPoints*maxEvalPointDist_ )
    {
        range = minNumEvalPoints*maxEvalPointDist_;
    }

    // find required number of evaluation points and corresponding step:
    size_t numEvalPoints = calculateNumEvalPoints(range);
    real deltaEvalPoints = range/(numEvalPoints - 1);

    // create set of evalPoints:
    real evalPointsLo = 0.5*(rangeHi + rangeLo) - 0.5*range;
    std::vector<real> evalPoints(numEvalPoints, evalPointsLo);
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
    numPoints = std::pow(2, std::ceil(std::log2(numPoints)));

    // return number of points:
    return numPoints;
}


/*!
 * Auxiliary function that carries out the actual kernel density estimation.
 * This is currently implemented as individual summations at each evaluation
 * point, i.e.
 *
 * \f[
 *      p(x) = \frac{1}{h N} \sum_{i=1}^{N} K\left( \frac{x - x_i}{h} \right)
 * \f]
 *
 * where \f$ h \f$ is the bandwidth, \f$ N \f$ is the number of samples, and
 * \f$ K(x) \f$ is a kernel function implemented as a class derived from
 * AbstractKernelFunction.
 *
 * Note that this is relatively costly for large sample sizes or many 
 * evaluation points and should in the future be replaced by a more efficient
 * FFT-based algorithm.
 *
 * \todo Implement the convolution as FFT. Note that Gromacs provides a wrapper
 * for Fourier transforms that should remain valid even if the library drops
 * FFTW.
 */
std::vector<real>
KernelDensityEstimator::calculateDensity(
        const std::vector<real> &samples,
        const std::vector<real> &evalPoints)
{
    // initialise density vector as zero:
    std::vector<real> density(evalPoints.size(), 0.0);

    // handle special case of empty sample:
    if( samples.size() == 0 )
    {
        // just return zero density:
        return density;
    }

    // create kernel:
    KernelFunctionPointer Kernel = KernelFunctionFactory::create(
            kernelFunction_);

    // normalisation constant:
    real normalisation = 1.0 / (samples.size() * bandWidth_);
    normalisation *= Kernel -> normalisingFactor();

    // scaled bandwidth:
    real bw = bandWidth_ * bandWidthScale_;

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints.size(); i++)
    {
        // density is sum over kernel distances:
        for(auto sample : samples)
        {
            density[i] += Kernel -> operator()( 
                    (evalPoints[i] - sample)/bw );
        }

        // normalise density at this evaluation point:
        density[i] *= normalisation;
    }

    // return density:
    return(density);
}


/*!
 * Auxiliary function for setting the density at the endpoints of the 
 * evaluation range to zero. This is done so that the SplineCurve1D returned
 * by evaluate() will return zero when extrapolating linearly. This in turn 
 * guarantees that the resulting spline curve will be integrable/normalisable. 
 *
 * Note that unless the kernel used for density estimation has only local 
 * support, the density will never be exactly zero, but the probability mass
 * far from the data range is generally negligible. The parameter set with
 * setEvalRangeCutoff() can be used to achieve a smooth transition between the
 * actual density and zero at the endpoints, as it controls how for the 
 * evaluation range extends beyond the data range (in multiples of the band
 * width). For normally distributed data, a cutoff of about three will ensure
 * that the probability mass neglected at the endpoints is less then 1%. Tests
 * indicate that a value of five is save for a wide variety of band widths and
 * evaluation steps.
 */
void
KernelDensityEstimator::endpointDensityToZero(
        std::vector<real> &density,
        std::vector<real> &evalPoints)
{
    // get evaluation point step length:
    real step = (evalPoints.back() - evalPoints.front())/(evalPoints.size() - 1);

    // append extra evaluation point at either side of the range:
    evalPoints.push_back(evalPoints.back() + step);
    evalPoints.insert(evalPoints.begin(), evalPoints.front() - step);

    // add zero density at either end of density vector:
    density.push_back(0.0);
    density.insert(density.begin(), 0.0);
}

