#ifndef NAIVE_CYLINDRICAL_PATH_FINDER_HPP
#define NAIVE_CYLINDRICAL_PATH_FINDER_HPP

#include <gromacs/utility/real.h> 
#include <gromacs/math/vec.h>

#include "path-finding/abstract_path_finder.hpp"


/*
 *
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

