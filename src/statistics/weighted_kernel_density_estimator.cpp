#include "geometry/linear_spline_interp_1D.hpp"
#include "statistics/weighted_kernel_density_estimator.hpp"



/*
 *
 */
SplineCurve1D
WeightedKernelDensityEstimator::estimate(
        std::vector<real> &samples,
        std::vector<real> &weights)
{    
    // sanity check:
    if( samples.size() != weights.size() )
    {
        throw std::logic_error("Need same number of samples and weights!");
    }

    // create evaluation points:
    std::vector<real> evalPoints = createEvaluationPoints(
            samples);

    // determine weighted density:
    std::vector<real> weightedDensity = calculateWeightedDensity(
            samples,
            weights,
            evalPoints);

    // set weighted endpoint density to zero:
    endpointDensityToZero(
            weightedDensity,
            evalPoints);

    // return weighted density as a spline curve:
    LinearSplineInterp1D Interp;
    return Interp(evalPoints, weightedDensity);
}


/*
 *
 */
std::vector<real>
WeightedKernelDensityEstimator::calculateWeightedDensity(
        std::vector<real> &samples,
        std::vector<real> &weights,
        std::vector<real> &evalPoints)
{
    // set up the density kerne:
    KernelFunctionPointer kernel = KernelFunctionFactory::create(
            kernelFunction_);

    // allocate the density vector:
    std::vector<real> density(evalPoints.size(), 0.0);
    std::vector<real> weightedDensity(evalPoints.size(), 0.0);


    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints.size(); i++)
    {
        // density is sum over kernel distances:
        for(size_t j = 0; j < samples.size(); j++)
        {
            // evaluate kernel function:
            real kern =  kernel -> operator()( 
                    evalPoints[i] - samples[j]/bandWidth_ );

            // for weighted and unweighted sums:
            density[i] += kern;
            weightedDensity[i] += kern*weights[j];
        }

        // Nadaraya-Watson estimate of local function value:
        weightedDensity[i] /= density[i];
    }

    // return density:
    return(weightedDensity);
}

