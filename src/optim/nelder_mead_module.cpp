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


#include <algorithm>

#include "optim/nelder_mead_module.hpp"


/*!
 * Constructor. Creates a NelderMeadModule object, but does not set any of its
 * properties.
 */
NelderMeadModule::NelderMeadModule()
{

}


/*!
 * Destructor.
 */
NelderMeadModule::~NelderMeadModule()
{

}


/*!
 * \brief Setter function for parameters.
 *
 * Takes a standard map as input, from which it tries to extract the parameters
 * required by the Nelder-Mead algorithm. Parameters are specified by using 
 * their name (as a string) as map key and their value as map value. 
 * Unrecognised entries will be ignored and an error will be thrown if required 
 * parameters without default are missing. Available options are:
 *
 *   - nmMaxIter: the maximum number of iterations to perform (required)
 *   - nmInitShift: shift of initial vertex coordinates with respect to guess point (required)
 *   - nmContractionPar: factor used in contraction step (defaults to 0.5)
 *   - nmExpansionPar: factor used in expansion step (defaults to 2.0)
 *   - nmReflectionPar: factor used in reflection step (defaults to 1.0)
 *   - nmShrinkagePar: factor used in shrinkage step (defaults to  0.5)
 */
void
NelderMeadModule::setParams(std::map<std::string, real> params)
{
    // number of iterations:
    if( params.find("nmMaxIter") != params.end() )
    {
        maxIter_ = params["nmMaxIter"];
    }
    else
    {
        std::cerr<<"ERROR: Maximum number of Nelder-Mead iterations not specified!"<<std::endl;
        std::abort();
    }

    // shift factor:
    if( params.find("nmInitShift") != params.end() )
    {
        initShiftFac_ = params["nmInitShift"];
    }
    else
    {
        std::cerr<<"ERROR: Shift factor for initial vertex generation not specified!"<<std::endl;
        std::abort();
    }

    // contraction parameter:
    if( params.find("nmContractionPar") != params.end() )
    {
        contractionPar_ = params["nmContractionPar"];
    }
    else
    {
        contractionPar_ = 0.5;
    }

    // expansion parameter:
    if( params.find("nmExpansionPar") != params.end() )
    {
        expansionPar_ = params["nmExpansionPar"];
    }
    else
    {
        expansionPar_ = 2.0;
    }

    // reflection parameter:
    if( params.find("nmReflectionPar") != params.end() )
    {
        reflectionPar_ = params["nmReflectionPar"];
    }
    else
    {
        reflectionPar_ = 1.0;
    }

    // reflection parameter:
    if( params.find("nmShrinkagePar") != params.end() )
    {
        shrinkagePar_ = params["nmShrinkagePar"];
    }
    else
    {
        shrinkagePar_ = 0.5;
    }
}


/*!
 * Sets the objective function to be maximised. An objective function takes a 
 * vector of reals as its only argument and returns a single real.
 */
void
NelderMeadModule::setObjFun(ObjectiveFunction objFun)
{
    this -> objFun_ = objFun;
}


/*!
 * Creates the initial simplex from one guess point. The first vertex of the 
 * initial simplex will simply be the guess point itself. All other vertices
 * are calculated by perturbing one coordinate of the guess vector, i.e.
 *
 * \f[
 *      \mathbf{x}_i = \mathbf{x}_1 + h\mathbf{e}_i
 * \f]
 *
 * where \f$\mathbf{e}_i\f$ is the \f$i\f$-th unit vector and \f$h\f$ is a 
 * small shift that can be set as a parameter. The objective function is not
 * evaluated at any vertex.
 */
void
NelderMeadModule::setInitGuess(std::vector<real> guess)
{
    // add guess point to simplex:
    OptimSpacePoint guessPoint;
    guessPoint.first = guess;
    simplex_.push_back(guessPoint);

    // add additional vertices to create a simplex of proper dimensionality:
    for(size_t i = 0; i < guess.size(); i++)
    {
        // copy guess point and perturb j-th component:
        OptimSpacePoint vertex = guessPoint;        
        vertex.first[i] += 1.0*initShiftFac_;

        // add vertex to simplex:
        simplex_.push_back(vertex);
    }
}


/*!
 * Performs the Nelder-Mead optimisation loop. Should only be called once
 * parameters, objective function, and initial point have been set.
 */
void
NelderMeadModule::optimise()
{
    // sanity checks:
    if( simplex_.size() != simplex_.back().first.size() + 1 )
    {
        std::cerr<<"ERROR: Incorrect number of simplex vertices!"<<std::endl;
        std::cerr<<"Did you call setInitGuess()?"<<std::endl;
        std::abort();
    }

    // initialise centroid as vector of all zeros:
    centroid_.first.insert(centroid_.first.begin(), simplex_.front().first.size(), 0.0); 

    // evaluate objective function at all vertices:
    std::vector<OptimSpacePoint>::iterator vert;
    for(vert = simplex_.begin(); vert != simplex_.end(); vert++)
    {
        vert -> second = objFun_(vert -> first);
    }

    // Nelder-Mead main loop:
    for(int i = 0; i < maxIter_; i++)
    {       
        // sort vertices by function values:
        std::sort(simplex_.begin(), 
                  simplex_.end(), 
                  CompOptimSpacePoints());
       
       
        // recalculate centroid:
        calcCentroid();

        // calculate the reflected point:
        OptimSpacePoint reflectedPoint = centroid_;
        reflectedPoint.scale(1.0 + reflectionPar_);
        reflectedPoint.addScaled(simplex_.front(), -reflectionPar_);

        // evaluate objective function at reflected point:
        reflectedPoint.second = objFun_(reflectedPoint.first);

        // reflected point better than second worst?
        if( comparison_(simplex_[1], reflectedPoint) )
        {
            // reflected point better than best?
            if( comparison_(simplex_.back(), reflectedPoint) )
            {
                // calculate expansion point:
                OptimSpacePoint expandedPoint = centroid_;
                expandedPoint.scale(1.0 - expansionPar_);
                expandedPoint.addScaled(reflectedPoint, expansionPar_);

                // evaluate objective function at expansion point:
                expandedPoint.second = objFun_(expandedPoint.first);

                // expanded point better than reflected point:
                if( expandedPoint.second < reflectedPoint.second )
                {
                    // accept expanded point:
                    simplex_.front() = expandedPoint;
                }
                else
                {
                    // accept reflected point:
                    simplex_.front() = reflectedPoint;
                }
            }
            else
            {
                // accept reflected point:
                simplex_.front() = reflectedPoint;
           }
        }
        else
        {
            // calculate contraction point: 
            OptimSpacePoint contractedPoint = centroid_;
            contractedPoint.scale(1.0 - contractionPar_);
            contractedPoint.addScaled(simplex_.front(), contractionPar_);

            // evaluate objective function at contracted point:
            contractedPoint.second = objFun_(contractedPoint.first);

            // contracted point better than worst?
            if( comparison_(simplex_.front(), contractedPoint) )
            {
                // accept contracted point:
                simplex_.front() = contractedPoint;
            }
            else
            { 
                // update all but the best vertex:
                std::vector<OptimSpacePoint>::iterator it;
                for(it = simplex_.begin(); it != simplex_.end() - 1; it++)
                {
                    // calculate shrinkage point:
                    it -> scale(shrinkagePar_);
                    it -> addScaled(simplex_.back(), 1.0 - shrinkagePar_);

                    // evaluate objective function at new vertex:
                    it -> second = objFun_(it -> first);
                }
            }
        }

        // ensure vertices are sorted:
        std::sort(simplex_.begin(), 
                  simplex_.end(), 
                  CompOptimSpacePoints());
    }
}


/*!
 * Returns the best point in optimisation space. Only meaningful if called
 * after optimise().
 */
OptimSpacePoint
NelderMeadModule::getOptimPoint()
{
    return simplex_.back();
}


/*!
 * Calculates the centroid of all except the first vertex of the simplex. This 
 * is usually the worst vertex.
 */
void
NelderMeadModule::calcCentroid()
{
        // reset centroid to vector of all zeros:
        std::fill(centroid_.first.begin(), centroid_.first.end(), 0.0);

        // loop over vertices and add to centroid_:
        std::vector<OptimSpacePoint>::iterator it;
        for(it = simplex_.begin() + 1; it != simplex_.end(); it++)
        {
            centroid_.add(*it);
        }

        // scale all elements by number:
        real fac = 1.0/(simplex_.size() - 1);
        centroid_.scale(fac);
}

