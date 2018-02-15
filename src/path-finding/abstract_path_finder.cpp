// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


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
 * Returns neighbourhood search cutoff.
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
 * Returns flag indicating whether neighbourhood search cutoff is set.
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
 * points and radii. Basically a wrapper around the constructor of 
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
AbstractPathFinder::setParameters(const PathFindingParameters& /*&params*/)
{
    parametersSet_ = true;
}

