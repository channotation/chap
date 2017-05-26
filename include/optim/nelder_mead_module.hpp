#ifndef NELDER_MEAD_MODULE_HPP
#define NELDER_MEAD_MODULE_HPP

#include <map>
#include <string>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>

#include "optim/optimisation.hpp"


/*! 
 * \brief Derivative-free optimisation using the Nelder-Mead downhill simplex
 * method.
 *
 * This class provides functionality for derivative-free optimisation of a
 * real-valued function in \f$N\f$ dimensions. After creating a 
 * NelderMeadModule, the algorithms parameters should be set with setParams() 
 * and an initial point needs to be specified with setInitGuess(). 
 * Additionally, the objective function must be set using setObjFun. After 
 * these setup steps, the optimisation procedure is started with optimise() and 
 * the best value found can be retrieved with getOptimPoint(). Note that this 
 * module performs maximisation rather than the canonical minimisation.
 *
 * The Nelder-Mead method is based on heuristically refining the vertex 
 * positions of a simplex in the optimisation space. In the \f$i\f$-th
 *
 *
 */
class NelderMeadModule
{
    friend class NelderMeadModuleTest;
    FRIEND_TEST(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest);

    public:

        // constructor and destructor:
        NelderMeadModule();
        ~NelderMeadModule();

        // setting parameters and initial point:
        void setParams(std::map<std::string, real> params);
        void setObjFun(ObjectiveFunction objFun);
        void setInitGuess(std::vector<real> guess);

        // optimisation and result retrieval:
        void optimise();
        OptimSpacePoint getOptimPoint();
        
    private:

        // control parameters:
        int maxIter_;
        real initShiftFac_;

        // internal parameters:
        real contractionPar_;
        real expansionPar_;
        real reflectionPar_;
        real shrinkagePar_;

        // objective function:
        ObjectiveFunction objFun_;

        // internal optimisation state:
        std::vector<OptimSpacePoint> simplex_;
        OptimSpacePoint centroid_; 

        // comparison functor:
        CompOptimSpacePoints comparison_;

        // internal functions:
        void calcCentroid();
};

#endif

