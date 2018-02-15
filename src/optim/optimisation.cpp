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


#include "optim/optimisation.hpp"


/******************************************************************************
 * OptimSpacePoint
 *****************************************************************************/

/*!
 * Adds another point to the point, but does not update the value of the
 * objective function at the new coordinates.
 */
void 
OptimSpacePoint::add(OptimSpacePoint other)
{
    for(size_t i = 0; i < this -> first.size(); i++)
    {
        this -> first[i] += other.first[i];
    }   
}


/*!
 * Adds another point to this point, where each coordinate of the other point
 * is scaled by a common factor. Does not update the function value at the new 
 * coordinate.
 */
void
OptimSpacePoint::addScaled(OptimSpacePoint other, real fac)
{
    for(size_t i = 0; i < this -> first.size(); i++)
    {
        this -> first[i] += other.first[i] * fac;
    }   
}


/*!
 * Scales the coordinates of the point by a given factor, but does not update 
 * the value of the objective function at the new coordinates.
 */
void
OptimSpacePoint::scale(real fac)
{
    for(size_t i = 0; i < this -> first.size(); i++)
    {
        this -> first[i] *= fac;
    }       
}


/*!
 * Returns the squared distance between the this and another optimisation space
 * point. Does not alter the internal state of either point.
 */
real
OptimSpacePoint::dist2(OptimSpacePoint other)
{
    real d = 0.0;
    for(size_t i = 0; i < this -> first.size(); i++)
    {
        d += (this->first[i] - other.first[i]) * (this->first[i] - other.first[i]);
    }
    return d;
}


/******************************************************************************
 * OptimisationModule
 *****************************************************************************/

/*!
 * Constructor.
 */
OptimisationModule::OptimisationModule()
{

}


/*!
 * Destructor.
 */
OptimisationModule::~OptimisationModule()
{

}

