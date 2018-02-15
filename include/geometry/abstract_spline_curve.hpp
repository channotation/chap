// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


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


/*!
 * \brief Abstract base class for spline curves in unspecified dimensions.
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
        std::vector<real> knots_;

        // basis spline (derivative) functor:
        BSplineBasisSet B_;

        // internal utility functions:
        int findInterval(const real &evalPoint);
};


#endif

