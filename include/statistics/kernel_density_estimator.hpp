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
            KernelDensityEstimatorParameterTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorEvalPointTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorGaussianRawDensityTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorGaussianInterpDensityTest);

    public:
        
        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // implementation of parameter setting interface:
        virtual void setParameters(
                const DensityEstimationParameters &params);

    private:

        // flag indicating whether parameters have been set:
        bool paramtersSet_;

        // internal parameters:
        real bandWidth_;
        real maxEvalPointDist_;
        real evalRangeCutoff_;
        eKernelFunction kernelFunction_;

        // auxiliary functions for parameter setting:
        void setBandWidth(const real bandWidth);
        void setMaxEvalPointDist(const real maxEvalPointDist);
        void setEvalRangeCutoff(const real evalRangeCutoff);
        void setKernelFunction(const eKernelFunction kernelFunction);

        // auxiliary functions for density estimation:
        std::vector<real> createEvaluationPoints(
                const std::vector<real> &samples);
        size_t calculateNumEvalPoints(
                const real range);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &evalPoints);
        void endpointDensityToZero(
                std::vector<real> &density,
                std::vector<real> &evalPoints);
};

#endif

