#ifndef KERNEL_DENSITY_ESTIMATOR_HPP
#define KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"
#include "statistics/abstract_density_estimator.hpp"
#include "statistics/kernel_function.hpp"


/*
 *
 */
class KernelDensityEstimator : public AbstractDensityEstimator
{
    friend class KernelDensityEstimatorTest;
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorEvalPointTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorGaussianDensityTest);

    public:
        
        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // implementation of parameter setting interface:
        virtual void setParameters(
                const DensityEstimationParameters &params);
        void setBandWidth(const real bandWidth);
        void setMaxEvalPointDist(const real maxEvalPointDist);
        void setEvalRangeCutoff(const real evalRangeCutoff);
        void setKernelFunction(const eKernelFunction kernelFunction);

    private:

        // internal parameters:
        real bandWidth_;
        real maxEvalPointDist_;
        real evalRangeCutoff_;
        eKernelFunction kernelFunction_;

        // auxiliary functions:
        std::vector<real> createEvaluationPoints(
                const std::vector<real> &samples);
        size_t calculateNumEvalPoints(
                const real range);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &evalPoints);
};

#endif

