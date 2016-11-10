#ifndef SIMULATEDANNEALINGMODULE_HPP
#define SIMULATEDANNEALINGMODULE_HPP

#include <gromacs/utility/real.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>


typedef	std::function<real(std::vector<real>)> energyFunction;


class simulatedAnnealingModule
{
	public:

	simulatedAnnealingModule(int stateDim,
							 int randomSeed,
							 int maxIter,
							 real initTemp,
							 real coolingFactor);
	
	void anneal();

	// getter functions (used in unit tests):
	int getStateDim(){return stateDim_;};
	int getMaxIter(){return maxIter_;};
	int getSeed(){return seed_;};
	real getTemp(){return temp_;};
	real getCoolingFactor(){return coolingFactor_;};

	private:

	int stateDim_;										// dimension of state space
	int maxIter_;										// maximum number of iterations
	int seed_;											// seed for random number generator
	
	real temp_;											// temperature
	real coolingFactor_;								// temperature reduction factor

	real crntEnergy_;									// cost function value of current state
	real candEnergy_;									// cost function value of candidate state
	real bestEnergy_;									// cost function value of best state

	std::vector<real> crntState_;						// state space coordinates of current state
	std::vector<real> candState_;						// state space coordinates of candidate state
	std::vector<real> bestState_;						// state space coordinates of best state

	gmx::DefaultRandomEngine rng_;						// pseudo random number generator
	gmx::UniformRealDistribution<real> candGenDistr_; 	// distribution for candidate generation
	gmx::UniformRealDistribution<real> candAccDistr_; 	// distribution for candidate acceptance

	energyFunction evaluateEnergy_;						// function object for evaluating energy of state

	// cool temperature:
	void cool();

	// generate a candidate state:
	void generateCandidate();

	// calculate probability of accepting candidate state:
	real calculateAcceptanceProbability();


};



#endif

