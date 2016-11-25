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
	FRIEND_TEST(SimulatedAnnealingModuleTest, IsotropicCandidateGenerationTest);
	FRIEND_TEST(SimulatedAnnealingModuleTest, AdaptiveCandidateGenerationTest);


	public:

	SimulatedAnnealingModule(int stateDim,
							 int randomSeed,
							 int maxCoolingIter,
							 int numCovSamples,
							 real initTemp,
							 real coolingFactor,
							 real stepLengthFactor,
							 real *initState,
							 costFunction cf,
							 bool useAdaptiveCandidateGeneration);

	~SimulatedAnnealingModule();

	void anneal();

	// getter functions (used in unit tests):
	int getStateDim(){return stateDim_;};
	int getMaxCoolingIter(){return maxCoolingIter_;};
	int getNumCovSamples(){return numCovSamples_;};
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

	bool useAdaptiveCandidateGeneration_;

	int stateDim_;										// dimension of state space
	int seed_;											// seed for random number generator
	int maxCoolingIter_;								// maximum number of cooling steps
	int numCovSamples_;									// candidate states generated per covariance matrix update

	int nCoolingIter_;									
	
	real temp_;											// temperature
	real coolingFactor_;								// temperature reduction factor
	real stepLengthFactor_;								// factor influencing candidate generation step

	real *crntState_;									// current state vector in optimisation space
	real *candState_;									// candidate state vector in optimisation space
	real *bestState_;									// best state vector in optimisation space

	real *stateSampleMatrix_;							// matrix containing state space sample collected in a cooling step
	real *covarianceMatrix_;							// covariance matrix
	real *adaptationMatrix_;							// adaptation matrix for candidate generation

	real crntCost_;										// cost function value at current state
	real candCost_;										// cost function value at candidate state
	real bestCost_;										// cost function value at best state

	gmx::DefaultRandomEngine rng_;						// pseudo random number generator
	gmx::UniformRealDistribution<real> candGenDistr_; 	// distribution for candidate generation
	gmx::UniformRealDistribution<real> candAccDistr_; 	// distribution for candidate acceptance

	costFunction evaluateCost_;


	void cool();
	void generateCandidateState();
	bool acceptCandidateState();
	void updateAdaptationMatrix();

};



#endif

