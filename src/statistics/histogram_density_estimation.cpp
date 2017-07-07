#include <algorithm>
#include <stdexcept>
#include <string>

#include "geometry/cubic_spline_interp_1D.hpp"
#include "geometry/linear_spline_interp_1D.hpp"

#include "statistics/histogram_density_estimator.hpp"

// TODO rmove
#include <iostream>

/*
 *
 */
HistogramDensityEstimator::HistogramDensityEstimator()
    : binWidth_(0.0)
{

}


/*
 *
 */
SplineCurve1D
HistogramDensityEstimator::estimate(
        std::vector<real> &samples)
{
    // sanity checks:
    if( binWidth_ <= 0 )
    {
        throw std::logic_error("Histogram bin width must be a positive number!");
    }

    // make sure input data is sorted:
    std::sort(samples.begin(), samples.end());

    // set up break points for this data set: 
    std::vector<real> breaks = createBreaks(
            samples.front(),
            samples.back());

    // compute midpoints corresponding to these breakpoints:
    std::vector<real> midpoints = createMidpoints(breaks);

    // emergency break for histogram size:
    // (If number of bins gets too large, the spline interpolation code will
    // fail with a segfault. This hardcoded limit is supposed to prevent that 
    // situation, but in practice the reasonable size limit for number of 
    // histogram bins is probably lower due to performance constraints.)
    int maxBinNumber = 25000;
    if( midpoints.size() > maxBinNumber )
    {
        throw std::runtime_error("Number of bins exceeds limit for spline "
        "interpolation! Need to increase bin width.");
    }
    
    // calculate density:
    std::vector<real> density = calculateDensity(samples, breaks);

    // sanity check:
    if( density.size() != midpoints.size() )
    {
        throw std::logic_error("Histogram has " + 
        std::string(std::to_string(density.size())) + " and " + 
        std::string(std::to_string(midpoints.size())) + " midpoints!");
    }

    // scale by inverse bin width to get proper density:
    std::transform(
            density.begin(),
            density.end(),
            density.begin(),
            std::bind2nd(std::multiplies<real>(), 1.0/binWidth_));

    // finally create spline curve from this:
    LinearSplineInterp1D Interp;
    return Interp(midpoints, density);
}


/*!
 * Setter method for assigning a bin width for the histogram.
 */
void
HistogramDensityEstimator::setBinWidth(
        real binWidth)
{
    // sanity check:
    if( binWidth <= 0.0 )
    {
        throw std::logic_error("Histogram bin width must be positive!");
    }

    // set internal bin width parameter:
    binWidth_ = binWidth;
}


/*
 *
 */
std::vector<real>
HistogramDensityEstimator::createBreaks(
        real rangeLo,
        real rangeHi)
{
    // will shift by half a bin width past lower endpoint:
    real halfBinWidth = 0.5*binWidth_;
    real breaksLo = rangeLo - 3.0*halfBinWidth;

    // build vector of breaks:
    std::vector<real> breaks;
    breaks.push_back(breaksLo);
    while( breaks.back() <= rangeHi + 3.0*halfBinWidth )
    {
        breaks.push_back(breaks.back() + binWidth_);
    }

    // return break points:
    return breaks;
}


/*
 *
 */
std::vector<real>
HistogramDensityEstimator::createMidpoints(
        const std::vector<real> &breaks)
{
    // reserve memory for midpoints:
    std::vector<real> midpoints;
    midpoints.resize(breaks.size() - 1);

    // midpoints are shifted wrt breaks by a half bin width:    
    std::transform(
            breaks.begin(),
            breaks.end() - 1,
            midpoints.begin(),
            std::bind2nd(std::plus<real>(), binWidth_));

    // return midpoints:
    return midpoints;
}


/*
 *
 */
std::vector<real>
HistogramDensityEstimator::calculateDensity(
        const std::vector<real> &samples,
        const std::vector<real> &breaks)
{
    // loop over intervals and count samples:
    std::vector<int> counts;
    counts.reserve(breaks.size() - 1);
    auto boundLo = samples.begin();
    auto boundHi = samples.end();
    for(size_t i = 0; i < breaks.size() - 1; i++)
    {
        // iterator to first element greater equal lower break point:
        boundLo = std::upper_bound(boundLo, samples.end(), breaks[i]);

        // iterator to first element strictly greater than upper break point:
        boundHi = std::upper_bound(boundLo, samples.end(), breaks[i+1]);

        // count number of samples in this range:
        counts.push_back(std::distance(boundLo, boundHi));
    }

    // sum over counts:
    int sum = std::accumulate(counts.begin(), counts.end(), 0.0);

    // sanity check:
    if( sum != samples.size() )
    {
        throw std::logic_error("Histogram counts sum to " +
        std::string(std::to_string(sum)) + " but number of samples is " +
        std::string(std::to_string(samples.size())) + "!");
    }

    // endpoint bins should be empty by construction:
    if( counts.front() != 0 || counts.back() != 0 )
    {
        throw std::logic_error("Histogram endpoint bins are not empty!");
    }

    // calculate density:
    std::vector<real> density;
    density.resize(counts.size());
    std::transform(
            counts.begin(),
            counts.end(),
            density.begin(),
            std::bind2nd(std::divides<real>(), sum));

    // return vector of densities:
    return(density);
}

