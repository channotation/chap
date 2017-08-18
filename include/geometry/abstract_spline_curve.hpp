#ifndef ABSTRACT_SPLINE_CURVE_HPP
#define ABSTRACT_SPLINE_CURVE_HPP

#include <vector>

#include <gromacs/math/vec.h>

#include "geometry/bspline_basis_set.hpp"


/*!
 * Shorthand notation for index zero in curvilinear coordinates for better
 * legibility.
 */
const short int SS = 0;

/*!
 * Shorthand notation for index zero in curvilinear coordinates for better
 * legibility.
 */
const short int RR = 1;

/*!
 * Shorthand notation for index zero in curvilinear coordinates for better
 * legibility.
 */
const short int PP = 2;


/*
 *
 */
class AbstractSplineCurve
{
    public:

        // getter methods:
        int degree() const;
        int nCtrlPoints() const;
        int nKnots() const;
        std::vector<real> knotVector() const;
        std::vector<real> uniqueKnots() const;

        // method to shift the internal coordinate system:
        void shift(const gmx::RVec &shift);


    protected:
        
        // internal variables:
        int degree_;
        int nCtrlPoints_;
        int nKnots_;
        std::vector<real> knotVector_;
        std::vector<real> knots_;

        // basis spline (derivative) functor:
        BSplineBasisSet B_;

        // internal utility functions:
        int findInterval(const real &evalPoint);
};


#endif

