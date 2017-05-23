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


/*
 * Setter function for parameters.
 */
void
NelderMeadModule::setParams(std::map<std::string, real> params)
{

}


/*
 *
 */
void
NelderMeadModule::optimise()
{
    maxIter_ = 100;


    for(int i = 0; i < maxIter_; i++)
    {
        // ORDERING
       
        // sort vertices by function values:
        std::sort(crntSimplex_.begin(), 
                  crntSimplex_.end(), 
                  CompOptimSpacePoints());
        
        

        // CENTROID


        // initialise centroid as vector of all zeros and infinite cost:
        OptimSpacePoint centroid;
        centroid.first.insert(centroid.first.begin(), crntSimplex_.front().first.size(), 0.0); 
        centroid.second = std::numeric_limits<real>::infinity();


        // loop over vertices and add to centroid:
        std::vector<OptimSpacePoint>::iterator it;
        for(it = crntSimplex_.begin(); it != crntSimplex_.end() - 1; it++)
        {
            centroid.add(*it);
        }

        // scale all elements by number:
        real fac = 1.0/(crntSimplex_.size() - 1);
        centroid.scale(fac);
        



        // REFLECTION

        


        // EXPANSION



        // CONTRACTION



        // SHRINK CONTRACTION




        // TERMINATION
        //
        //
    }
}

