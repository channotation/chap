#ifndef KERNEL_DENSITY_ESTIMATOR_HPP
#define KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"

#include "statistics/abstract_density_estimator.hpp"
#include "statistics/kernel_function.hpp"


/*
 *
 */
class KernelDensityEstimator : public AbstractDensityEstimator
{
    public:
        
        // density estimation interface:
        virtual SplineCurve1D estimate(
                const std::vector<real> &samples);

        // implementation of parameter setting interface:
        virtual void setParameters(
                const DensityEstimationParameters &params);
        void setBandWidth(const real bandWidth);
        void setMaxEvalPointDist(const real maxEvalPointDist);
        void setKernelFunction(const eKernelFunction kernelFunction);

    private:

        // internal parameters:
        real bandWidth_;
        real maxEvalPointDist_;
        eKernelFunction kernelFunction_;



        //
        std::vector<real> createEvaluationPoints(
                const std::vector<real> &samples);
        size_t calculateNumEvalPoints(
                const real range);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &evalPoints);
};

#endif

