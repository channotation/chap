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


#ifndef SPLINE_CURVE_1D_HPP
#define SPLINE_CURVE_1D_HPP

#include <utility>
#include <vector>

#include <gtest/gtest_prod.h>   

#include <gromacs/math/vec.h>

#include "geometry/abstract_spline_curve.hpp"


/*!
 * \brief Spline curve in one dimension.
 *
 * This class represents a spline curve in one spatial dimension, i.e. a spline
 * function. In three dimensions, the class SplineCurve3D can be used.
 */
class SplineCurve1D : public AbstractSplineCurve
{
    friend class SplineCurve1DTest;
    FRIEND_TEST(SplineCurve1DTest, SplineCurve1DFindIntervalTest);

    public:
    
        // constructor and destructor:
        SplineCurve1D(
                int degree, 
                std::vector<real> knotVector,
                std::vector<real> ctrlPoints); 
        SplineCurve1D();

        // public interface for curve evaluation:
        real evaluate(
                const real &eval, 
                unsigned int deriv);
        std::vector<real> evaluateMultiple(
                const std::vector<real> &eval, 
                unsigned int deriv);

        // getter function for control points:
        std::vector<real> ctrlPoints() const;

        // compute spline properties:
        real length() const;
        std::pair<real, real> minimum(const std::pair<real, real> &lim);

    private:

        // internal variables:
        std::vector<real> ctrlPoints_;

        // auxiliary functions for evaluation:
        inline real evaluateInternal(const real &eval, unsigned int deriv);
        inline real evaluateExternal(const real &eval, unsigned int deriv);
        inline real computeLinearCombination(const SparseBasis &basis);
};

#endif

