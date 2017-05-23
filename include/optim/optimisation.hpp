#ifndef OPTIMISATION_HPP
#define OPTIMISATION_HPP

#include <utility>
#include <vector>

#include <gromacs/utility/real.h>


/*!
 * Representation of a point in optimisation space as a pair of a vector of
 * coordinates in the optimisation space and the corresponding value of the
 * objective function at this point.
 */
class OptimSpacePoint : public std::pair<std::vector<real>, real>
{
    public:

        void add(OptimSpacePoint other);
        void scale(real fac);
};


/*!
 * Functor for comparing two points in optimisation space by their respective
 * function value. 
 */
typedef struct CompOptimSpacePoints
{
    bool operator()(OptimSpacePoint pointA, OptimSpacePoint pointB)
    {
        return pointB.second < pointB.second;
    }
} CompOptimSpacePoints;

#endif

