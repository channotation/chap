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


#ifndef ABSTRACT_PATH_FINDER_HPP
#define ABSTRACT_PATH_FINDER_HPP

#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/molecular_path.hpp"


/*!
 * \brief Enaum for available path-finding methods.
 */
typedef enum {ePathFindingMethodNaiveCylindrical,
              ePathFindingMethodInplaneOptimised} ePathFindingMethod;


/*!
 * \brief Helper class for specifying parameters in the classes derived from
 * AbstractPathFinder.
 *
 * This class is used to simplify the interface of AbstractPathFinder by 
 * providing a simple method of passing parameters to this class. All getter
 * methods perform a check for whether the requested parameter has been set.
 */
class PathFindingParameters
{
    public:

        // constructor and destructor:
        PathFindingParameters();

        // setter methods:
        void setNbhCutoff(real nbhCutoff);
        void setProbeStepLength(real probeStepLength);
        void setMaxProbeRadius(real maxProbeRadius);
        void setMaxProbeSteps(int maxProbeSteps);

        // getter methods:
        real nbhCutoff() const;
        bool nbhCutoffIsSet() const;
        
        real probeStepLength() const;
        bool probeStepLengthIsSet() const;

        real maxProbeRadius() const;
        bool maxProbeRadiusIsSet() const;

        int maxProbeSteps() const;
        bool maxProbeStepsIsSet() const;

    private:

        real nbhCutoff_;
        bool nbhCutoffIsSet_;

        real probeStepLength_;
        bool probeStepLengthIsSet_;

        real maxProbeRadius_;
        bool maxProbeRadiusIsSet_;

        int maxProbeSteps_;
        bool maxProbeStepsIsSet_;
};


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

        // public interface for path retrieval:
        virtual MolecularPath getMolecularPath();

        // public interface for setting parameters:
        virtual void setParameters(const PathFindingParameters& /* &params*/);

        // convenience functions for retrieving path points and radii directly:
        std::vector<gmx::RVec> pathPoints(){return path_;};
        std::vector<real> pathRadii(){return radii_;};


    protected:

        // internal map of parameters:
        std::map<std::string, real> params_;

        // have parameters been set:
        bool parametersSet_;
    
        // data containers for path points and corresponding radii:
        std::vector<gmx::RVec> path_;
        std::vector<real> radii_;
};

#endif

