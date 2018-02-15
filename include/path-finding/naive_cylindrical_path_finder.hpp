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

