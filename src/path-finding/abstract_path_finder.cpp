#include <map>
#include <string>

#include "path-finding/abstract_path_finder.hpp"


/*!
 * Constructor sets parameter values to nonsensical values and flags to false.
 */
PathFindingParameters::PathFindingParameters()
    : probeStepLength_(-1.0)
    , probeStepLengthIsSet_(false)
    , maxProbeRadius_(-1.0)
    , maxProbeRadiusIsSet_(false)
    , maxProbeSteps_(0)
    , maxProbeStepsIsSet_(false)
{

}


/*!
 * Sets cutoff for neighborhood search.
 */
void
PathFindingParameters::setNbhCutoff(real nbhCutoff)
{
    nbhCutoff_ = nbhCutoff;
    nbhCutoffIsSet_ = true;
}


/*!
 * Sets probe step length. 
 */
void
PathFindingParameters::setProbeStepLength(real probeStepLength)
{
    probeStepLength_ = probeStepLength;
    probeStepLengthIsSet_ = true;
}


/*!
 * Sets maximum probe radius.
 */
void
PathFindingParameters::setMaxProbeRadius(real maxProbeRadius)
{
    maxProbeRadius_ = maxProbeRadius;
    maxProbeRadiusIsSet_ = true;
}


/*!
 * Sets maximum number of probe steps.
 */
void
PathFindingParameters::setMaxProbeSteps(int maxProbeSteps)
{
    maxProbeSteps_ = maxProbeSteps;
    maxProbeStepsIsSet_ = true;
}


/*!
 * Returns neighborhodd search cutoff.
 *
 * \throws std::logic_error If parameter value unset/
 */
real
PathFindingParameters::nbhCutoff() const
{
    if( nbhCutoffIsSet_ )
    {
        return nbhCutoff_;
    }
    else
    {
        throw std::logic_error("Parameter nbhCutoff is not set.");
    }
}


/*!
 * Returns flag indicating whether neighborhodd search cutoff is set.
 */
bool
PathFindingParameters::nbhCutoffIsSet() const
{
    return nbhCutoffIsSet_;
}


/*!
 * Returns probe step length.
 *
 * \throws std::logic_error If parameter value unset/
 */
real
PathFindingParameters::probeStepLength() const
{
    if( probeStepLengthIsSet_ )
    {
        return probeStepLength_;
    }
    else
    {
        throw std::logic_error("Parameter probeStepLength is not set.");
    }
}


/*!
 * Returns flag indicating if probe step length has been set.
 */
bool
PathFindingParameters::probeStepLengthIsSet() const
{
    return probeStepLengthIsSet_;
}


/*!
 * Returns maximum probe radius.
 *
 * \throws std::logic_error If parameter value unset/
 */
real
PathFindingParameters::maxProbeRadius() const
{
    if( maxProbeRadiusIsSet_ )
    {
        return maxProbeRadius_;
    }
    else
    {
        throw std::logic_error("Parameter maxProbeRadius is not set.");
    }
}


/*!
 * Returns flag indicating whether maximum probe radius has been set.
 */
bool
PathFindingParameters::maxProbeRadiusIsSet() const
{
    return maxProbeRadiusIsSet_;
}


/*!
 * Returns maximum number of probe steps.
 *
 * \throws std::runtime_error If parameter value unset/ 
 */
int
PathFindingParameters::maxProbeSteps() const
{
    if( maxProbeStepsIsSet_ )
    {
        return maxProbeSteps_;
    }
    else
    {
        throw std::logic_error("Parameter maxProbeSteps is not set.");
    }
}


/*!
 * Returns flag indicating if maximum number of probe steps has been set.
 */
bool
PathFindingParameters::maxProbeStepsIsSet() const
{
    return maxProbeStepsIsSet_;
}



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


/*!
 * Setter function for parameters. On the base class, only changes the flag.
 */
void
AbstractPathFinder::setParameters(const PathFindingParameters &params)
{
    parametersSet_;
}


