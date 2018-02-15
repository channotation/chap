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


#ifndef CUBIC_SPLINE_INTERP_1D_HPP
#define CUBIC_SPLINE_INTERP_1D_HPP

#include <vector>

#include "geometry/abstract_cubic_spline_interp.hpp"
#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Functor class for performing cubic spline interpolation in one 
 * dimension.
 */
class CubicSplineInterp1D : public AbstractCubicSplineInterp
{
    public:
        
        // constructor and destructor:
        CubicSplineInterp1D();
        ~CubicSplineInterp1D();

        // interpolation interface:
        SplineCurve1D interpolate(std::vector<real> &x,
                                  std::vector<real> &f,
                                  eSplineInterpBoundaryCondition bc);
        SplineCurve1D operator()(std::vector<real> &x,
                                 std::vector<real> &f,
                                 eSplineInterpBoundaryCondition bc);

        // curve properties:
        std::pair<real, real> findMinimum() const;
};

#endif

