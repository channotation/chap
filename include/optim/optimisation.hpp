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


#ifndef OPTIMISATION_HPP
#define OPTIMISATION_HPP

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gromacs/utility/real.h>


/*!
 * \brief Representation of a point in optimisation space and its corresponding
 * objective function value.
 *
 * Representation of a point in optimisation space as a pair of a vector of
 * coordinates in the optimisation space and the corresponding value of the
 * objective function at this point. The class also provides a few convenience
 * functions for manipulating optimisation space point coordinates that are
 * utilised by the various optimisation modules.
 */
class OptimSpacePoint : public std::pair<std::vector<real>, real>
{
    public:

        void add(OptimSpacePoint other);
        void addScaled(OptimSpacePoint other, real fac);
        void scale(real fac);

        real dist2(OptimSpacePoint other);
};


/*!
 * \brief Functor for comparing two points in optimisation space by their
 * respective function value. 
 *
 * Returns true if the value of the objective function is lower at the first 
 * point than at the second point. This allows sorting a vector of
 * OptimSpacePoints using standard library functions.
 */
typedef struct CompOptimSpacePoints
{
    bool operator()(OptimSpacePoint pointA, OptimSpacePoint pointB)
    {
        return pointA.second < pointB.second;
    }
} CompOptimSpacePoints;


/*!
 * \typedef Shorthand notation for an objective function as used by 
 * the Nelder-Mead optimisation method.
 */
typedef std::function<real(std::vector<real>)> ObjectiveFunction;


/*!
 * \brief Abstract base class for optimisation modules.
 *
 * Specifies public interface that needs to be implemented by an optimisation 
 * module.
 */
class OptimisationModule
{
    public:

        // constructor and destructor:
        OptimisationModule();
        ~OptimisationModule();

        // public interface for optimisation classes:
        virtual void setParams(std::map<std::string, real>) = 0;
        virtual void setObjFun(ObjectiveFunction objFun) = 0;
        virtual void setInitGuess(std::vector<real> guess) = 0;
        virtual void optimise() = 0;
        virtual OptimSpacePoint getOptimPoint() = 0;
};

#endif

