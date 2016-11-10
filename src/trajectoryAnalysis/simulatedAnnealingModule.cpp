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
												   real coolingFactor)
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
{
	// initialise state vectors:
	crntState_ = {0.0, 0.0, 0.0}; // initialise properly! --> pass as argument!
	candState_ = crntState_;
	bestState_ = crntState_;

	// get initial energy of states:
	crntEnergy_ = evaluateEnergy_(crntState_);
	candEnergy_ = crntEnergy_;
	bestEnergy_ = bestEnergy_;
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
		
		// evaluate cost function
		candEnergy_ = evaluateEnergy_(candState_);

		// calculate acceptance probability:
		accProb = calculateAcceptanceProbability();

		// accept move?
		if( accProb < candAccDistr_(rng_) )
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
	return std::exp( (candEnergy_ - crntEnergy_)/temp_ );
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
 *
 */
void
simulatedAnnealingModule::generateCandidate()
{

}


