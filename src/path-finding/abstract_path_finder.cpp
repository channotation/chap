#include <map>
#include <string>

#include "path-finding/abstract_path_finder.hpp"


/*
 * Constructor to be used in initialiser list of derived classes.
 */
AbstractPathFinder::AbstractPathFinder(std::map<std::string, real> params)
    : path_()
    , radii_()
    , params_(params)
{

}


/*
 *
 */
MolecularPath
AbstractPathFinder::getMolecularPath()
{
    return MolecularPath(path_, radii_);
}

