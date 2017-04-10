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
    public:

  
    protected:

        // internal variables:
        int degree_;
        int nCtrlPoints_;
        int nKnots_;
        std::vector<real> knotVector_;

        // basis spline functor:
        BasisSpline B_;

        // internal utility functions:
        int findInterval(real &evalPoint);
        real deBoorRecursion(int r, 
                             int i, 
                             real &evalPoint, 
                             const std::vector<real> &ctrlCoefs);
};


#endif

