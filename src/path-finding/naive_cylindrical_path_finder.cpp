#include <cmath>

#include "path-finding/naive_cylindrical_path_finder.hpp"


/*
 * Constructor.
 */
NaiveCylindricalPathFinder::NaiveCylindricalPathFinder(real stepLength,
                                                       int nSteps,
                                                       real cylRad,
                                                       gmx::RVec centrePoint,
                                                       gmx::RVec dirVec)
{
    // normalise direction vector:
    gmx::RVec normDirVec = dirVec;
    real norm = std::sqrt( dirVec[0]*dirVec[0] + 
                           dirVec[1]*dirVec[1] +
                           dirVec[2]*dirVec[2] );
    normDirVec[0] = normDirVec[0]/norm;
    normDirVec[1] = normDirVec[1]/norm;
    normDirVec[2] = normDirVec[2]/norm;

    // build up path points and radii:
    for(int i = 0; i <= 2*nSteps; i++)
    {
        real x = centrePoint[0] + (i - nSteps)*stepLength*normDirVec[0];
        real y = centrePoint[1] + (i - nSteps)*stepLength*normDirVec[1];
        real z = centrePoint[2] + (i - nSteps)*stepLength*normDirVec[2];
        real r = cylRad;

        path_.push_back(gmx::RVec(x, y, z));
        radii_.push_back(cylRad);
    }
}


/*
 * Destructor.
 */
NaiveCylindricalPathFinder::~NaiveCylindricalPathFinder()
{

}


/*
 *
 */
void
NaiveCylindricalPathFinder::findPath()
{

}

