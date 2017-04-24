#ifndef SPLINE_CURVE_3D_HPP
#define SPLINE_CURVE_3D_HPP

#include <vector>

#include <gtest/gtest_prod.h>   

#include <gromacs/utility/real.h> 
#include <gromacs/math/vec.h>

#include "geometry/abstract_spline_curve.hpp"


/*
 *
 */
class SplineCurve3D : public AbstractSplineCurve
{
    public:
        
        // constructor and destructor:
        SplineCurve3D(int degree,
                      std::vector<real> knotVector,
                      std::vector<gmx::RVec> ctrlPoints);
        ~SplineCurve3D();

        // public interface for spline evaluation:
        gmx::RVec evaluate(real &evalPoint,
                           unsigned int derivOrder,
                           eSplineEvalMethod method);
        gmx::RVec operator()(real &evalPoint,
                             unsigned int derivOrder,
                             eSplineEvalMethod method);

        // reparameterisation methods:
        void arcLengthParam();

        // calculate properties of curve:
        real length(real &lo, real &hi, const real &absTol);
        real length(const real &absTol);
 
        gmx::RVec tangentVec(real &evalPoint);
        

        real speed(real &evalPoint);
        
        
    private:

        // internal variables:
        std::vector<gmx::RVec> ctrlPoints_;

        // utility functions:
        real arcLengthBoole(real &lo, real &hi);


};

#endif
