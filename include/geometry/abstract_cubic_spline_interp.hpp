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


#ifndef ABSTRACT_CUBIC_SPLINE_INTERP_HPP
#define ABSTRACT_CUBIC_SPLINE_INTERP_HPP

#include <vector>

#include <gromacs/math/vec.h>

#include "geometry/basis_spline.hpp"


enum eSplineInterpBoundaryCondition {eSplineInterpBoundaryHermite, 
                                     eSplineInterpBoundaryNatural};    
enum eSplineInterpEndpoint {eSplineInterpEndpointLo, 
                            eSplineInterpEndpointHi};
enum eSplineInterpDerivEstimate {eSplineInterpDerivSimple,
                                 eSplineInterpDerivParabolic};


/*!
 * \brief Abstract base class providing common functionality for cubic spline 
 * interpolation in one and three dimensions.
 *
 * This class provides the methods used internally by CubicSplineInterp1D and 
 * CubicSplineInterp3D. In particular it implements the methods used for 
 * assembling the tridiagonal system matrix and right hand side vector defining
 * a functional interpolation problem. In addition, it provides methods for
 * generating an appropriate knot vector from a set of abscissa points and a 
 * simple estimate of the input function's derivative that is needed for some 
 * types of boundary conditions.
 *
 * Input data is assumed to be given in the form of \f$ N \f$ abscissa points 
 * \f$ x_i \f$ and corresponding function values \f$ f_i = f(x_i) \f$. From 
 * this, a cubic spline \f$ \sigma(x) \f$ that smoothly interpolates between all 
 * input data points can be found by solving a tridiagonal system of equations:
 *
 * \f[
 *      \mathbf{Ac} =
 *      \begin{bmatrix}
 *          \alpha_{1} & \gamma_{1} &             &             &               \\
 *          \beta_{2}  & \alpha_{2} & \gamma_{2}  &             &               \\
 *                     & \ddots     & \ddots      & \ddots      &               \\
 *                     &            & \beta_{N+1} & \alpha_{N+1} & \gamma_{N+1} \\
 *                     &            &             & \beta_{N+2}  & \alpha_{N+2} 
 *      \end{bmatrix}
 *      \begin{bmatrix}
 *          c_{1} \\
 *          c_{2} \\
 *          \vdots \\
 *          c_{N + 1} \\
 *          c_{N + 2} 
 *      \end{bmatrix}
 *      =
 *      \begin{bmatrix}
 *          f'(x_1) \\
 *          f(x_1) \\
 *          \vdots \\
 *          f(x_N) \\
 *          f'(x_N) 
 *      \end{bmatrix}
 *      = \mathbf{f}
 * \f]
 *
 * Here \f$ \mathbf{c} \f$ is the vector of control points of 
 * \f$ \sigma(x) \f$ and the corresponding knot vector is given by
 *
 * \f[ 
 *      \{t_i\}_{i=1}^{N+6} = \{x_1, x_1, x_1, x_1, x_2, ..., x_{N-1}, x_N, x_N, x_N, x_N\} 
 * \f]
 *
 * The second and second-to-last rows in the above system represent the 
 * Hermite boundary conditions \f$ \sigma'(x_1) = f'(x_1) \f$ and 
 * \f$ \sigma'(x_N) = f'(x_N) \f$, while the intermediate rows represent the
 * interpolation condition \f$ \sigma(x_i) = f(x_i) \f$. Consequently, the 
 * nonzero elements of the system matrix are given by
 *
 * \f[
 *      \alpha_i = B_{i, 3}(x_{i-1}) \\
 *      \beta_i = B_{i-1, 3}(x_{i-1}) \\
 *      \gamma_i = B_{i+1, 3}(x_{i-1})
 * \f]
 *
 * for \f$ 2 \leq i \leq N+1 \f$ and 
 *
 * \f[
 *      \alpha_1 = B'_{1, 3}(x_1) \\
 *      \alpha_{N+2} = B'_{N+2, 3}(x_{N}) \\
 *      \beta_{N+2} = B'_{N+1, 3}(x_{N}) \\
 *      \gamma_1 = B'_{2, 3}(x_1) 
 * \f]
 *
 * at the boundaries. The evaluation of the 
 * relevant basis splines and their derivatives is handled by the BasisSpline
 * and BasisSplineDerivative functors respectively. The derivatives of
 * \f$ f(x) \f$ occurring in the right hand side vector can be approximated by a
 * simple finite difference
 *
 * \f[
 *      f'(x_1) \approx f(x_2) - f(x_1) \quad \text{and} \\
 *      f'(x_N) \approx f(x_N) - f(x_{N-1})
 * \f]
 *
 * or alternatively by a parabolic approximation of \f$ f(x) \f$,
 *
 * \f[
 *      f'(x_1) \approx \frac{\Delta x_0 \delta_1 + \Delta x_1 \delta_0}{\Delta x_0 + \Delta x_1} \quad \text{and} \\
 *      f'(x_N) \approx \frac{\Delta x_{N-1} \delta_{N} + \Delta x_{N} \delta_{N-1}}{\Delta x_{N-1} + \Delta x_{N}}
 * \f]
 *
 * can be used. Here \f$ \Delta x_i = x_{i+1} - x_i \f$ and 
 *
 * \f[
 *      \delta_i = \frac{f(x_{i+1}) - f(x_i))}{\Delta x_i}
 * \f]
 *
 * and the ghost points \f$ x_0 = x_3 \f$ and \f$ x_{N+1} = x_{N-1} \f$ are 
 * introduced.
 *
 * AbstractCubicSplineInterp provides the utilities for correctly assembling 
 * the system matrix and right hand side, the routines for solving the system 
 * are implemented in the derived classes CubicSplineInterp1D and 
 * CubicSplineInterp3D.
 */
class AbstractCubicSplineInterp
{
    public:

        // constructor and destructor:
        AbstractCubicSplineInterp();
        ~AbstractCubicSplineInterp();

    protected:

        // member variables:
        const int degree_ = 3;
        eSplineInterpBoundaryCondition bc_;

        // internal helper functions:
        void assembleDiagonals(std::vector<real> &knotVector,
                               std::vector<real> &x,
                               real *subDiag,
                               real *mainDiag,
                               real *superDiag,
                               eSplineInterpBoundaryCondition bc);
        void assembleRhs(std::vector<real> &x,
                         std::vector<real> &f,
                         real *rhsVec,
                         eSplineInterpBoundaryCondition bc);
        std::vector<real> prepareKnotVector(std::vector<real> &x);
        real estimateEndpointDeriv(std::vector<real> &x,
                                   std::vector<real> &f,
                                   eSplineInterpEndpoint endpoint,
                                   eSplineInterpDerivEstimate method);
};

#endif

