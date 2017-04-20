#ifndef ABSTRACT_SPLINE_CURVE_HPP
#define ABSTRACT_SPLINE_CURVE_HPP

#include <vector>

#include <gromacs/math/vec.h>

#include "geometry/basis_spline.hpp"


/*
 * Enum for spline evaluation method.
 */
enum eSplineEvalMethod {eSplineEvalNaive = 901, eSplineEvalDeBoor = 902};


/*
 *
 */
class AbstractSplineCurve
{ 
    protected:

        // internal variables:
        int degree_;
        int nCtrlPoints_;
        int nKnots_;
        std::vector<real> knotVector_;

        // basis spline (derivative) functor:
        BasisSpline B_;
        BasisSplineDerivative D_;

        // internal utility functions:
        int findInterval(real &evalPoint);
        real deBoorRecursion(int r, 
                             int i, 
                             real &evalPoint, 
                             const std::vector<real> &ctrlCoefs);

        // internal drivers for evaluation method:
        real evaluateSplineFun(real &evalPoint,
                               const std::vector<real> &ctrlCoefs,
                               unsigned int derivOrder, 
                               eSplineEvalMethod method);
        real evaluateNaive(real &evalPoint,
                           const std::vector<real> &ctrlCoefs);
        real evaluateDeBoor(real &evalPoint,
                            const std::vector<real> &ctrlCoefs);
        real evaluateDeriv(real &evalPoint,
                           const std::vector<real> &ctrlCoefs, 
                           unsigned int order);
        
        // internal drivers for extrapolation method:
        real linearExtrap(real &evalPoint,
                          const std::vector<real> ctrlPoints,
                          real &boundary, 
                          unsigned int derivOrder);
};


#endif

