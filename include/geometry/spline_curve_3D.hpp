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
        SplineCurve3D();
        SplineCurve3D(const SplineCurve3D &other);
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
        real length(real &lo, real &hi);
        real length();
     
        gmx::RVec tangentVec(real &evalPoint);
//        gmx::RVec normalVec(real &evalPoints);
        
        real speed(real &evalPoint);

        // utilities for accessing arc length at the control points:
        std::vector<real> ctrlPointArcLength();
        real frstPointArcLength();
        real lastPointArcLength();

        // map points onto curve:
        real pointSqDist(gmx::RVec point, real eval);
        int closestCtrlPoint(gmx::RVec &point);
        gmx::RVec cartesianToCurvilinear(gmx::RVec cartPoint,
                                         real lo,
                                         real hi,
                                         real tol);
 
        // getter functions:
        std::vector<gmx::RVec> ctrlPoints() const;
        
    private:

        // internal variables:
        std::vector<gmx::RVec> ctrlPoints_;

        bool arcLengthTableAvailable_;
        std::vector<real> arcLengthTable_;

        // curve length utilities:
        real arcLengthBoole(real lo, real hi);
        void prepareArcLengthTable();
        
        // arc length reparameterisation utilities:
        real arcLengthToParam(real &arcLength);
        bool arcLengthToParamTerm(real lo, real hi, real tol);
        real arcLengthToParamObj(real lo, real hi, real target);
};

#endif

