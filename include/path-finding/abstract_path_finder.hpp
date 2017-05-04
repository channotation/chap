#ifndef ABSTRACT_PATH_FINDER_HPP
#define ABSTRACT_PATH_FINDER_HPP

#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/molecular_path.hpp"


/*
 * Abstract base class for all permeation path finding algorithms. Should 
 * provide the public interface.
 */
class AbstractPathFinder
{
    public:

        // constructor:
        AbstractPathFinder();
       
        // interface for path finding method:
        virtual void findPath() = 0;

        // public interface for path retrieal:
        MolecularPath getMolecularPath();


//        virtual std::vector<gmx::RVec> getPath(){return path_;};
//        virtual std::vector<real> getRadii(){return radii_;};


    protected:

 
        std::vector<gmx::RVec> path_;
        std::vector<real> radii_;

};

#endif

