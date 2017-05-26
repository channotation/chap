#include <algorithm>

#include "optim/nelder_mead_module.hpp"







/*
 * Constructor.
 */
NelderMeadModule::NelderMeadModule()
{

}


/*
 * Destructor.
 */
NelderMeadModule::~NelderMeadModule()
{

}


/*!
 * \brief Setter function for parameters.
 *
 * Takes a standard map as input, from which it tries to extract the parameters
 * required by the Nelder-Mead algorithm. Unrecognised entries will be ignored
 * and an error will be thrown if required parameters without default are 
 * missing.
 */
void
NelderMeadModule::setParams(std::map<std::string, real> params)
{
    // number of iterations:
    if( params.find("maxIter") != params.end() )
    {
        maxIter_ = params["maxIter"];
    }
    else
    {
        std::cerr<<"ERROR: Maximum number of Nelder-Mead iterations not specified!"<<std::endl;
        std::abort();
    }

    // shift factor:
    if( params.find("initShift") != params.end() )
    {
        initShiftFac_ = params["initShift"];
    }
    else
    {
        std::cerr<<"ERROR: Shift factor for initial vertex generation not specified!"<<std::endl;
        std::abort();
    }

    // contraction parameter:
    if( params.find("contractionPar") != params.end() )
    {
        contractionPar_ = params["contractionPar"];
    }
    else
    {
        contractionPar_ = 0.5;
    }

    // expansion parameter:
    if( params.find("expansionPar") != params.end() )
    {
        expansionPar_ = params["expansionPar"];
    }
    else
    {
        expansionPar_ = 2.0;
    }

    // reflection parameter:
    if( params.find("reflectionPar") != params.end() )
    {
        reflectionPar_ = params["reflectionPar"];
    }
    else
    {
        reflectionPar_ = 1.0;
    }

    // reflection parameter:
    if( params.find("shrinkagePar") != params.end() )
    {
        shrinkagePar_ = params["shrinkagePar"];
    }
    else
    {
        shrinkagePar_ = 0.5;
    }
}


/*
 *
 */
void
NelderMeadModule::setObjFun(ObjectiveFunction objFun)
{
    this -> objFun_ = objFun;
}


/*
 *
 */
void
NelderMeadModule::setInitGuess(std::vector<real> guess)
{
    // add guess point to simplex:
    OptimSpacePoint guessPoint;
    guessPoint.first = guess;
    simplex_.push_back(guessPoint);

    // add additional vertices to create a simplex of proper dimensionality:
    for(int i = 0; i < guess.size(); i++)
    {
        // copy guess point and perturb j-th component:
        OptimSpacePoint vertex = guessPoint;        
        vertex.first[i] += 1.0*initShiftFac_;

        // add vertex to simplex:
        simplex_.push_back(vertex);
    }
}


/*
 *
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

    // initialse centroid as vector of all zeros:
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
 * Returns the best point in optimisation space. Only meaningfull if called
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

