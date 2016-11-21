#include <iostream>

#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"


/*
 * Constructor.
 *
 * Note that the candidate generation range is [-sqrt(3), sqrt(3)), which gives
 * a standard deviation of 1.0 for a uniform distribution.
 */
simulatedAnnealingModule::simulatedAnnealingModule(int stateDim,
                                                   int randomSeed,
												   int maxIter,
												   real initTemp,
												   real coolingFactor,
												   std::vector<real> initState,
												   energyFunction ef)
	: stateDim_(stateDim)
	, seed_(randomSeed)
	, maxIter_(maxIter)
	, temp_(initTemp)
	, coolingFactor_(coolingFactor)
	, rng_(seed_)
	, candGenDistr_(-1.732051, 1.732051)
	, candAccDistr_(0.0, 1.0)
	, crntState_()
	, candState_()
	, bestState_()
	, evaluateEnergy_(ef)
{
	// initialise state vectors:
	crntState_ = initState;
	candState_ = initState;
	bestState_ = initState;

	// get initial energy of states:
	crntEnergy_ = evaluateEnergy_(crntState_);
	candEnergy_ = evaluateEnergy_(candState_);
	bestEnergy_ = evaluateEnergy_(bestState_);
}


/*
 * Perform annealing procedure.
 *
 * TODO: Handle constraints!
 */
void
simulatedAnnealingModule::anneal()
{
	// initialise temporary variables:
	real accProb = 0.0;

	// simulated annealing iteration:
	for(int i = 0; i < maxIter_; i++)
	{

		// generate candidate
		generateCandidate();

		// evaluate cost function
		candEnergy_ = evaluateEnergy_(candState_);

		// calculate acceptance probability:
		accProb = calculateAcceptanceProbability();
/*
		// inform user:
		std::cout<<"i = "<<i
				 <<"  T = "<<temp_
				 <<"  P = "<<accProb
		         <<"  crntState = "<<crntState_[0]<<" , "<<crntState_[1]
				 <<"  crntEnergy = "<<crntEnergy_
		         <<"  candState = "<<candState_[0]<<" , "<<candState_[1]
				 <<"  candEnergy = "<<candEnergy_
				 <<"  bestState = "<<bestState_[0]<<" , "<<bestState_[1]
				 <<"  bestEnergy = "<<bestEnergy_<<std::endl;
*/
		// accept move?
		if( candAccDistr_(rng_) < accProb )
		{
			// candidate state becomes current state:
			crntEnergy_ = candEnergy_;
			crntState_ = candState_;

			// new best state found?
			if( candEnergy_ > bestEnergy_ )
			{
				// candidate state becomes best state:
				bestEnergy_ = candEnergy_;
				bestState_ = candState_;
			}
		}

		// reduce temperature:
		cool();
	}
}


/*
 * Calculate acceptance probability from Boltzmann distribution.
 * TODO: Should cost be capped at 1?
 */
real
simulatedAnnealingModule::calculateAcceptanceProbability()
{
	real cap = 1.0;
	return std::min(std::exp( (candEnergy_ - crntEnergy_)/temp_ ), cap);
}



/*
 * Reduces temperature of SA module. Currently only simple exponential 
 * cooling is implemented.
 *
 * TODO: Implement additional cooling schedules.
 */
void
simulatedAnnealingModule::cool()
{
	temp_ *= coolingFactor_;
}



/*
 * Generates a candidate state from the current state.
 */
void
simulatedAnnealingModule::generateCandidate()
{
	// loop over state dimension:
	for(int i=0; i<stateDim_; i++)
	{
		// new state is current state plus some small random offset:
		// TODO: introduce step length:
		candState_.at(i) = crntState_.at(i) + 0.01*candGenDistr_(rng_);
	}
}


