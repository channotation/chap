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
    for(int i = 0; i < this -> first.size(); i++)
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
    for(int i = 0; i < this -> first.size(); i++)
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
    for(int i = 0; i < this -> first.size(); i++)
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
    for(int i = 0; i < this -> first.size(); i++)
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
