#ifndef OPTIMISATION_HPP
#define OPTIMISATION_HPP

#include <functional>
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
 * Functor for comparing two points in optimisation space by their respective
 * function value. Returns true if the value of the objective function is lower
 * at the first point than at the second point. This allows sorting a vector of
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
 * Shorthand notation for an objective function as used by the Nelder-Mead
 * optimisation method.
 */
typedef std::function<real(std::vector<real>)> ObjectiveFunction;


#endif

