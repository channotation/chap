#ifndef SIMULATEDANNEALINGMODULE_HPP
#define SIMULATEDANNEALINGMODULE_HPP

#include <gtest/gtest_prod.h>

#include <gromacs/utility/real.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>


typedef std::function<real(real*)> costFunction;


class SimulatedAnnealingModule
{
	friend class SimulatedAnnealingModuleTest;
	FRIEND_TEST(SimulatedAnnealingModuleTest, CoolingTest);


	public:

	SimulatedAnnealingModule(int stateDim,
							 int randomSeed,
							 int maxIter,
							 real initTemp,
							 real coolingFactor,
							 real stepLengthFactor,
							 real *initState,
							 costFunction cf);

	~SimulatedAnnealingModule();

	void anneal();

	// getter functions (used in unit tests):
	int getStateDim(){return stateDim_;};
	int getMaxIter(){return maxIter_;};
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

	bool useAdaptiveCandidateGeneration;

	int stateDim_;										// dimension of state space
	int maxIter_;										// maximum number of iterations
	int seed_;											// seed for random number generator
	int maxCoolingIter_;								// maximum number of cooling steps
	int candGenTrials_;									// candidate states generated per covariance matrix update

	int nCoolingIter_;
	
	real temp_;											// temperature
	real coolingFactor_;								// temperature reduction factor
	real stepLengthFactor_;								// factor influencing candidate generation step

	real *crntState_;
	real *candState_;
	real *bestState_;

	real crntCost_;
	real candCost_;
	real bestCost_;



//	real crntEnergy_;									// cost function value of current state
//	real candEnergy_;									// cost function value of candidate state
//	real bestEnergy_;									// cost function value of best state

//	std::vector<real> crntState_;						// state space coordinates of current state
//	std::vector<real> candState_;						// state space coordinates of candidate state
//	std::vector<real> bestState_;						// state space coordinates of best state

	gmx::DefaultRandomEngine rng_;						// pseudo random number generator
	gmx::UniformRealDistribution<real> candGenDistr_; 	// distribution for candidate generation
	gmx::UniformRealDistribution<real> candAccDistr_; 	// distribution for candidate acceptance

//	energyFunction evaluateEnergy_;						// function object for evaluating energy of state
	costFunction evaluateCost_;

//	void generateStateDirection();
//	bool acceptCandidate();

	// cool temperature:
	void cool();

	void generateCandidateState();

	// generate a candidate state:
//	void generateCandidate();

	// calculate probability of accepting candidate state:
//	real calculateAcceptanceProbability();
};



#endif

