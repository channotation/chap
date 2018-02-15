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


#ifndef NELDER_MEAD_MODULE_HPP
#define NELDER_MEAD_MODULE_HPP

#include <map>
#include <string>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>

#include "optim/optimisation.hpp"


/*! 
 * \brief Derivative-free optimisation using the Nelder-Mead downhill simplex
 * method.
 *
 * This class provides functionality for derivative-free optimisation of a
 * real-valued function in \f$ D \f$ dimensions. After creating a 
 * NelderMeadModule, the algorithms parameters should be set with setParams() 
 * and an initial point needs to be specified with setInitGuess(). 
 * Additionally, the objective function must be set using setObjFun(). After 
 * these setup steps, the optimisation procedure is started with optimise() and 
 * the best value found can be retrieved with getOptimPoint(). Note that this 
 * module performs maximisation rather than the canonical minimisation.
 *
 * The Nelder-Mead method is based on heuristically refining the vertex 
 * positions of a simplex in the optimisation space. A simplex is in 
 * \f$ D \f$-dimensional space is a geometric figure defined by \f$ D+1 \f$
 * vertices, i.e. a triangle in two dimensions and a tetrahedron in three
 * dimensions. Starting from some initial guess, the vertices are moved through
 * optimisation space until the simplex finally shrinks around the position of
 * the optimum (which may be local rather than global).
 *
 * In each iteration of the Nelder-Mead algorithm, all vertex positions are 
 * sorted by their respective objective function value. Subsequently, four
 * different operations, namely reflection, expansion, contraction, and 
 * shrink contraction may be carried out to adjust the vertex positions. 
 * Letting \f$ \mathbf{v}_w \f$, \f$ \mathbf{v}_s \f$, and \f$ \mathbf{v}_b \f$
 * denote the worst, second-worst, and best vertex respectively, a reflection
 * point is generated as
 *
 * \f[
 *    \mathbf{v}_r = \mathbf{c} + \alpha \left( \mathbf{c} - \mathbf{v}_w \right)
 * \f]
 *
 * where \f$ \mathbf{c} \f$ denotes the centre of geometry of all but the worst
 * vertex and \f$ \alpha \f$ is a parameter that is usually set to 1. If the
 * objective function has a higher value at the reflection point than at the 
 * second worst current vertex, but not higher than the current best vertex
 * (i.e. \f$ f_s < f_r < f_b \f$ ), than the reflection point replaces the 
 * current worst vertex and the next iteration of the algorithm begins.
 *
 * If the reflection point improves even on the current best vertex 
 * \f$ f_r > f_b \f$ than the algorithm tries to explore this direction in 
 * optimisation space even further and an expansion point is calculated 
 * according to 
 *
 * \f[
 *    \mathbf{v}_e = \mathbf{c} + \gamma \left( \mathbf{v}_r - \mathbf{c} \right)
 * \f]
 *
 * where the parameter \f$ \gamma \f$ is commonly set to 2. The better point
 * out of \f$ \mathbf{v}_r \f$ and \f$ \mathbf{v}_e \f$ will then replace
 * \f$ \mathbf{v}_w \f$ and the next iteration begins.
 *
 * However, if the reflection point is worse than the current second worst
 * vertex (\f$ f_r < f_s \f$ ), a contraction point is calculated as
 *
 * \f[
 *    \mathbf{v}_c = \mathbf{c} + \beta \left( \mathbf{v}_w - \mathbf{c} \right)
 * \f]
 *
 * with the parameter \f$ \beta \f$ typically set to 0.5. If the contraction 
 * point is an improvement (\f$ f_c > f_w \f$ ), it will be accepted as new 
 * vertex. 
 *
 * Otherwise a shrink contraction is carried out and all but the best current
 * vertex are recalculated according to
 *
 * \f[
 *    \mathbf{v}_j + \mathbf{v}_b + \delta \left( \mathbf{v}_j - \mathbf{v}_b \right)
 * \f]
 *
 * where the parameter \f$ \delta \f$ is typically set to 0.5. The algorithm is
 * terminated after a maximum number of iterations.
 */
class NelderMeadModule : public OptimisationModule
{
    friend class NelderMeadModuleTest;
    FRIEND_TEST(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest);

    public:

        // constructor and destructor:
        NelderMeadModule();
        ~NelderMeadModule();

        // setting parameters and initial point:
        void setParams(std::map<std::string, real> params);
        void setObjFun(ObjectiveFunction objFun);
        void setInitGuess(std::vector<real> guess);

        // optimisation and result retrieval:
        void optimise();
        OptimSpacePoint getOptimPoint();
        
    private:

        // control parameters:
        int maxIter_;
        real initShiftFac_;

        // internal parameters:
        real contractionPar_;
        real expansionPar_;
        real reflectionPar_;
        real shrinkagePar_;

        // objective function:
        ObjectiveFunction objFun_;

        // internal optimisation state:
        std::vector<OptimSpacePoint> simplex_;
        OptimSpacePoint centroid_; 

        // comparison functor:
        CompOptimSpacePoints comparison_;

        // internal functions:
        void calcCentroid();
};

#endif

