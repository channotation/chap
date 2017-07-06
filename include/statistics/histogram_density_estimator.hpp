#ifndef HISTOGRAM_DENSITY_ESTIMATOR_HPP
#define HISTOGRAM_DENSITY_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "statistics/abstract_density_estimator.hpp"


/*
 *
 */
class HistogramDensityEstimator : public AbstractDensityEstimator
{
    friend class HistogramDensityEstimatorTest;
    FRIEND_TEST(HistogramDensityEstimatorTest, HistogramDensityEstimatorBreaksTest);
    FRIEND_TEST(HistogramDensityEstimatorTest, HistogramDensityEstimatorEstimateTest);

    public:
       
        // constructor and destructor:
        HistogramDensityEstimator();

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // setter methods:
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

};

#endif

