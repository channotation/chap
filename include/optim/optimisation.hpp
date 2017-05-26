#ifndef OPTIMISATION_HPP
#define OPTIMISATION_HPP

#include <functional>
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
        void addScaled(OptimSpacePoint other, real fac);
        void scale(real fac);

        real dist2(OptimSpacePoint other);
};


/*!
 * Functor for comparing two points in optimisation space by their respective
 * function value. 
 */
typedef struct CompOptimSpacePoints
{
    bool operator()(OptimSpacePoint pointA, OptimSpacePoint pointB)
    {
        return pointA.second < pointB.second;
    }
} CompOptimSpacePoints;


/*
 *
 */
typedef std::function<real(std::vector<real>)> ObjectiveFunction;


#endif

