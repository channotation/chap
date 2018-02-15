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


#ifndef NAIVE_CYLINDRICAL_PATH_FINDER_HPP
#define NAIVE_CYLINDRICAL_PATH_FINDER_HPP

#include <gromacs/utility/real.h> 
#include <gromacs/math/vec.h>

#include "path-finding/abstract_path_finder.hpp"


/*!
 * \brief Fallback path-finding module, which constructs a simple cylinder as
 * molecular pathway.
 *
 * This is not technically finding a path, but may be useful if other 
 * path-finding algorithms fail and the user can manually define a cylindrical
 * pathway.
 */
class NaiveCylindricalPathFinder : public AbstractPathFinder
{
    public:

        // constructor and destructor:
        NaiveCylindricalPathFinder(std::map<std::string, real> params, 
                                   gmx::RVec centrePoint,
                                   gmx::RVec dirVec);
        ~NaiveCylindricalPathFinder();

        // path finding interface:
        virtual void findPath();

    private:

        int nSteps_;
        real cylRad_;
        real stepLength_;
        gmx::RVec centrePoint_;
        gmx::RVec dirVec_;

};


#endif

