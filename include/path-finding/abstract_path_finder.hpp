#ifndef ABSTRACT_PATH_FINDER_HPP
#define ABSTRACT_PATH_FINDER_HPP

#include <vector>

#include <gromacs/trajectoryanalysis.h>

/*
 * Abstract base class for all permeation path finding algorithms. Should 
 * provide the public interface.
 */
class AbstractPathFinder
{
    public:
       
        virtual void findPath() = 0;

        std::vector<gmx::RVec> getPath(){return path_;};
        std::vector<real> getRadii(){return radii_;};


    protected:

        AbstractPathFinder();
 
        std::vector<gmx::RVec> path_;
        std::vector<real> radii_;

};

#endif

