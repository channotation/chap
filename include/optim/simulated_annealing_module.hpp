#ifndef SIMULATED_ANNEALING_MODULE_HPP
#define SIMULATED_ANNEALING_MODULE_HPP

#include <map>
#include <string>

#include <gtest/gtest_prod.h>

#include <gromacs/utility/real.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>

#include "optim/optimisation.hpp"


/*!
 * \brief Multidimensional optimisation using simulated annealing.
 *
 * This class implements a simple version of the classic simulated annealing
 * algorithm for multidimensional optimisation. 
 *
 * \todo Document parameters properly. 
 */
class SimulatedAnnealingModule : public OptimisationModule
{
	friend class SimulatedAnnealingModuleTest;

	public:

        // constructors and destructors:
        SimulatedAnnealingModule();
		~SimulatedAnnealingModule();

        // public interface:
        virtual void setParams(std::map<std::string, real> params);
        virtual void setObjFun(ObjectiveFunction objFun);
        virtual void setInitGuess(std::vector<real> objFun);
        virtual void optimise();
        OptimSpacePoint getOptimPoint();

		// getter functions (used in unit tests):
		int getStateDim(){return stateDim_;};
		int getMaxCoolingIter(){return maxCoolingIter_;};
		int getSeed(){return seed_;};

		real getTemp(){return temp_;};
		real getCoolingFactor(){return coolingFactor_;};
		real getStepLengthFactor(){return stepLengthFactor_;};

		real getCrntStateAt(int i){return crntState_[i];};
		real getCandStateAt(int i){return candState_[i];};
		real getBestStateAt(int i){return bestState_[i];};

		std::vector<real> getCrntState(){return crntState_;};
		std::vector<real> getCandState(){return candState_;};
		std::vector<real> getBestState(){return bestState_;};

		real getCrntCost(){return crntCost_;};
		real getCandCost(){return candCost_;};
		real getBestCost(){return bestCost_;};

	private:

        // internal driver for annealing procedure:
        void anneal();

		// parameters:
		int seed_;	    		// seed for random number generator
		int stateDim_;	    	// dimension of state space
		int maxCoolingIter_;	// maximum number of cooling steps
	
		// internal state variables:
		real temp_;		    	// temperature
		real coolingFactor_;	// temperature reduction factor
		real stepLengthFactor_;	// factor for candidate generation step

        std::vector<real> crntState_;
        std::vector<real> candState_;
        std::vector<real> bestState_;
		//real *crntState_;		// current state in optimisation space
		//real *candState_;		// candidate state in optimisation space
		//real *bestState_;		// best state vector in optimisation space

		real crntCost_; 		// cost function value at current state
		real candCost_;    	// cost function value at candidate state
		real bestCost_;			// cost function value at best state

        // random number generation:
		gmx::DefaultRandomEngine rng_;		
		gmx::UniformRealDistribution<real> candGenDistr_; 	
        gmx::UniformRealDistribution<real> candAccDistr_; 	

		// functors and function type members:
        ObjectiveFunction objFun_;

		// member functions
		void annealIsotropic();
		void cool();
		void generateCandidateStateIsotropic();
		bool acceptCandidateState();
};

#endif

