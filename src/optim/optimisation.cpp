#include "optim/optimisation.hpp"


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

