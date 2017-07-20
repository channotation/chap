#include <map>
#include <string>

#include "path-finding/abstract_path_finder.hpp"


/*!
 * \brief Constructor to be used in initialiser list of derived classes. 
 *
 * Initialises containers for path points and radii and sets internal map of
 * parameters.
 */
AbstractPathFinder::AbstractPathFinder(std::map<std::string, real> params)
    : params_(params)
    , path_()
    , radii_()
{

}


/*!
 * Returns a MolecularPath object constructed from the path finder's path 
 * points and radii. Basically jsut a wrapper around the constructor of 
 * MolecularPath.
 */
MolecularPath
AbstractPathFinder::getMolecularPath()
{
    return MolecularPath(path_, radii_);
}

