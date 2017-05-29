#ifndef SIMULATED_ANNEALING_MODULE_HPP
#define SIMULATED_ANNEALING_MODULE_HPP

#include <map>
#include <string>

#include <gtest/gtest_prod.h>

#include <gromacs/utility/real.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>

#include "optim/optimisation.hpp"
#include "statistics-utilities/calculate_covariance_matrix.hpp"


typedef std::function<real(real*)> costFunction;


enum eSimAnTerm {CONVERGENCE = 101, 
                 MAX_COOLING_ITER = 102, 
                 NO_CAND_ACCEPTED = 103};


class SimulatedAnnealingModule
{
	friend class SimulatedAnnealingModuleTest;
	FRIEND_TEST(SimulatedAnnealingModuleTest, CoolingTest);
	FRIEND_TEST(SimulatedAnnealingModuleTest, IsotropicCandidateGenerationTest);
	FRIEND_TEST(SimulatedAnnealingModuleTest, AdaptiveCandidateGenerationTest);

	public:

		SimulatedAnnealingModule(int stateDim,
								 int randomSeed,
								 int maxCoolingIter,
								 int numCovSamples,
								 real xi,
								 real convRelTol,
								 real initTemp,
								 real coolingFactor,
								 real stepLengthFactor,
								 real *initState,
								 costFunction cf,
								 bool useAdaptiveCandidateGeneration);

        SimulatedAnnealingModule();
		~SimulatedAnnealingModule();

		eSimAnTerm anneal();

        // public interface:
        void setParams(std::map<std::string, real> params);
        void setObjFun(ObjectiveFunction objFun);
        void setInitGuess(std::vector<real> objFun);
        void optimise();
        OptimSpacePoint getOptimPoint();

		// getter functions (used in unit tests):
		int getStateDim(){return stateDim_;};
		int getMaxCoolingIter(){return maxCoolingIter_;};
		int getNumCandSamples(){return numCostSamples_;};
		int getSeed(){return seed_;};

		real getTemp(){return temp_;};
		real getCoolingFactor(){return coolingFactor_;};
		real getStepLengthFactor(){return stepLengthFactor_;};

		real getCrntStateAt(int i){return crntState_[i];};
		real getCandStateAt(int i){return candState_[i];};
		real getBestStateAt(int i){return bestState_[i];};

		real* getCrntState(){return crntState_;};
		real* getCandState(){return candState_;};
		real* getBestState(){return bestState_;};

		real getCrntCost(){return crntCost_;};
		real getCandCost(){return candCost_;};
		real getBestCost(){return bestCost_;};

	private:

		// parameters:
		bool useAdaptiveCandidateGeneration_;

		int seed_;									// seed for random number generator
		int stateDim_;								// dimension of state space
		int maxCoolingIter_;							// maximum number of cooling steps
		int numCostSamples_;							// candidate states generate per check of convergence criterion
	
		real beta_;									// free random walk parameter from Vanderbilt & Louie
		real xi_;										// growth factor from Vanderbilt & Louie
		real convRelTol_;

		// internal state variables:
		real temp_;											// temperature
		real coolingFactor_;								// temperature reduction factor
		real stepLengthFactor_;								// factor influencing candidate generation step

		real *crntState_;									// current state vector in optimisation space
		real *candState_;									// candidate state vector in optimisation space
		real *bestState_;									// best state vector in optimisation space

		real *stateSampleMatrix_;							// matrix containing state space sample collected in a cooling step
		real *adaptationMatrix_;							// adaptation matrix for candidate generation

		real crntCost_;										// cost function value at current state
		real candCost_;										// cost function value at candidate state
		real bestCost_;										// cost function value at best state

		real *costSamples_;

		gmx::DefaultRandomEngine rng_;						// pseudo random number generator
		gmx::UniformRealDistribution<real> candGenDistr_; 	// distribution for candidate generation
		gmx::UniformRealDistribution<real> candAccDistr_; 	// distribution for candidate acceptance

		// funcotrs and function type members:
		costFunction evaluateCost;
        ObjectiveFunction objFun_;
		CalculateCovarianceMatrix calculateCovarianceMatrix;

		// member functions
		eSimAnTerm annealIsotropic();						// function for non-adaptive annealing with isotropic canidate generation
		eSimAnTerm annealAdaptive();						// function for adaptive annealing

		void cool();
		void generateCandidateStateIsotropic();
		void generateCandidateStateAdaptive();
		void updateAdaptationMatrix();

		bool acceptCandidateState();						// checks Boltzmann criterion for accepting new candidate
		bool isConvergedIsotropic();						// checks if the non-adaptive algorithm has reached convergence
		bool isConvergedAdaptive();							// checks if the adaptive algorithm has reached convergence
};



#endif

