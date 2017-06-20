#ifndef ABSTRACT_PATH_FINDER_HPP
#define ABSTRACT_PATH_FINDER_HPP

#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/molecular_path.hpp"


/*!
 * \brief Abstract base class for all permeation path finding algorithms. 
 *
 * This specifies an interface for all path finding classes, namely that they
 * should provide a method findPath that runs the path finding algorithm.
 */
class AbstractPathFinder
{
    public:

        // constructor:
        AbstractPathFinder(std::map<std::string, real> params);
       
        // interface for path finding method:
        virtual void findPath() = 0;

        // public interface for path retrieal:
        virtual MolecularPath getMolecularPath();

        // convenience functions for retrieving path points and radii directly:
        std::vector<gmx::RVec> pathPoints(){return path_;};
        std::vector<real> pathRadii(){return radii_;};


    protected:

        // internal map of parameters:
        std::map<std::string, real> params_;
 
        // data containers for path points and corresponding radii:
        std::vector<gmx::RVec> path_;
        std::vector<real> radii_;
};

#endif

