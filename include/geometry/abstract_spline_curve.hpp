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

#include "geometry/basis_spline.hpp"


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
 * Enum for spline evaluation method.
 */
enum eSplineEvalMethod {eSplineEvalNaive = 901, eSplineEvalDeBoor = 902};


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

//    protected:
        
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

