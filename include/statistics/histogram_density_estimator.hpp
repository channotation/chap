#ifndef HISTOGRAM_DENSITY_ESTIMATOR_HPP
#define HISTOGRAM_DENSITY_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "statistics/abstract_density_estimator.hpp"


/*!
 * \brief One-dimensional density estimate using a histogram.
 *
 * This class provides a method estimate() for estimating a probability 
 * density of from a sample of scalar values. It uses a histogram with
 * uniformly spaced bins, where the bin width can be set using setBinWidth().
 *
 * The resulting density is then interpolated linearly using 
 * LinearSplineInterp1D and the result is returned as a SplineCurve1D object.
 * Linear rather than higher order interpolation is used to avoid artifacts
 * where the interpolation method introduces locally negative densities.
 */
class HistogramDensityEstimator : public AbstractDensityEstimator
{
    friend class HistogramDensityEstimatorTest;
    FRIEND_TEST(HistogramDensityEstimatorTest, HistogramDensityEstimatorBreaksTest);
    FRIEND_TEST(HistogramDensityEstimatorTest, HistogramDensityEstimatorDensityTest);
    FRIEND_TEST(HistogramDensityEstimatorTest, HistogramDensityEstimatorEstimateTest);

    public:
       
        // constructor and destructor:
        HistogramDensityEstimator();

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // setter methods for parameters:
        virtual void setParameters(
                const DensityEstimationParameters &params);
        void setBinWidth(real binWidth);

    private:

        // internal parameters:
        real binWidth_;

        // auxiliary functions:
        std::vector<real> createBreaks(
                real rangeLo,
                real rangeHi);
        std::vector<real> createMidpoints(
                const std::vector<real> &breaks);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &breaks);
    
};

#endif

