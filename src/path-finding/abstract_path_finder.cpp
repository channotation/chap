#include "path-finding/abstract_path_finder.hpp"

/*
 * Constructor to be used in initialiser list of derived classes.
 */
AbstractPathFinder::AbstractPathFinder()
    : path_()
    , radii_()
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

