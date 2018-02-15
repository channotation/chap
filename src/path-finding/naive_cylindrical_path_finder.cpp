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


#include <cmath>
#include <iostream>

#include "path-finding/naive_cylindrical_path_finder.hpp"


/*!
 * Constructor.
 */
NaiveCylindricalPathFinder::NaiveCylindricalPathFinder(std::map<std::string, real> params,
                                                       gmx::RVec centrePoint,
                                                       gmx::RVec dirVec)
    : AbstractPathFinder(params)
    , centrePoint_(centrePoint)
    , dirVec_(dirVec)
{
    // assign parameters:
    if( params.find("pfCylStepLength") != params.end() )
    {
        stepLength_ = params["pfCylStepLength"];
    }
    else
    {
        std::cerr<<"ERROR: Step length not given!"<<std::endl;
        std::abort();
    }

    if( params.find("pfCylRad") != params.end() )
    {
        cylRad_ = params["pfCylRad"];
    }
    else
    {
        std::cerr<<"ERROR: Cylinder radius not given!"<<std::endl;
        std::abort();
    }

    if( params.find("pfCylNumSteps") != params.end() )
    {
        nSteps_ = params["pfCylNumSteps"];
    }
    else
    {
        std::cerr<<"ERROR: Number of steps not given!"<<std::endl;
        std::abort();
    }
}


/*!
 * Destructor.
 */
NaiveCylindricalPathFinder::~NaiveCylindricalPathFinder()
{

}


/*!
 * Build cylindrical pathway.
 */
void
NaiveCylindricalPathFinder::findPath()
{
    // normalise direction vector:
    gmx::RVec normDirVec = dirVec_;
    real norm = std::sqrt( dirVec_[XX]*dirVec_[XX] + 
                           dirVec_[YY]*dirVec_[YY] +
                           dirVec_[ZZ]*dirVec_[ZZ] );
    normDirVec[XX] = normDirVec[XX]/norm;
    normDirVec[YY] = normDirVec[YY]/norm;
    normDirVec[ZZ] = normDirVec[ZZ]/norm;

    // build up path points and radii:
    for(int i = 0; i <= 2*nSteps_; i++)
    {
        real x = centrePoint_[XX] + (i - nSteps_)*stepLength_*normDirVec[XX];
        real y = centrePoint_[YY] + (i - nSteps_)*stepLength_*normDirVec[YY];
        real z = centrePoint_[ZZ] + (i - nSteps_)*stepLength_*normDirVec[ZZ];

        path_.push_back(gmx::RVec(x, y, z));
        radii_.push_back(cylRad_);
    }
}

