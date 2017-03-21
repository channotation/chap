#ifndef SIMULATED_ANNEALING_MODULE_HPP
#define SIMULATED_ANNEALING_MODULE_HPP

#include <gtest/gtest_prod.h>

#include <gromacs/utility/real.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>

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

		~SimulatedAnnealingModule();

		eSimAnTerm anneal();

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
		const bool useAdaptiveCandidateGeneration_;

		const int seed_;									// seed for random number generator
		const int stateDim_;								// dimension of state space
		const int maxCoolingIter_;							// maximum number of cooling steps
		const int numCostSamples_;							// candidate states generate per check of convergence criterion
	
		const real beta_;									// free random walk parameter from Vanderbilt & Louie
		const real xi_;										// growth factor from Vanderbilt & Louie
		const real convRelTol_;

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

