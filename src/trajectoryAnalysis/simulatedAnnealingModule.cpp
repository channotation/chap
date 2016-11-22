#include <iostream>

#include <cblas.h>

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
												   real stepLengthFactor,
												   std::vector<real> initState,
												   energyFunction ef)
	: stateDim_(stateDim)
	, seed_(randomSeed)
	, maxIter_(maxIter)
	, temp_(initTemp)
	, coolingFactor_(coolingFactor)
	, stepLengthFactor_(stepLengthFactor)
	, rng_(seed_)
	, candGenDistr_(-1.732051, 1.732051)
	, candAccDistr_(0.0, 1.0)
	, crntState_()
	, candState_()
	, bestState_()
//	, evaluateEnergy_(ef)
{
	// initialise state vectors:
//	crntState_ = new crnt
//	candState_ = initState;
//	bestState_ = initState;

	// get initial energy of states:
//	crntEnergy_ = evaluateEnergy_(crntState_);
//	candEnergy_ = evaluateEnergy_(candState_);
//	bestEnergy_ = evaluateEnergy_(bestState_);

	// get cost of inital states:
//	crntCost_ = evaluateCost_(crntState_);
///	candCost_ = evaluateCost_(candState_);
//	bestCost_ = evaluateCost_(bestState_);
}


/*
 * Perform annealing procedure.
 *
 * TODO: Handle constraints!
 */
void
simulatedAnnealingModule::anneal()
{
/*
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
*/
	///////////////////////////////////////////////////////////////////////////
	
	// start annealing loop:
	while(true)
	{
		// reset counter of accepted and rejected moves:
		int nAccepted = 0;
		int nRejected = 0;

		// 
		for(int i = 0; i < candGenTrials_; i++)
		{
			// generate a candidate state:
			generateCandidateState();

			// TODO: check boundary conditions!
			
			// evaluate cost function:
			candCost_ = evaluateCost_(candState_);

			// accept candidate?
			if( true )
			{
				// candidate state becomes current state:
				cblas_scopy(stateDim_, crntState_, 1, candState_, 1);
				crntCost_ = candCost_;

				// increment acceptance counter:
				nAccepted++;

				// is new state also the best state?
				if( crntCost_ > bestCost_ )
				{
					cblas_scopy(stateDim_, bestState_, 1, candState_, 1);
					bestCost_ = candCost_;
				}

			}
			else
			{
				// increment rejection counter:
				nRejected++;
			}

		}
	

		// maximum step number reached?
		if( nCoolingIter_ >= maxCoolingIter_ )
		{
			break;
		}

		// no candidates accepted?
		if( nAccepted <= 0 )
		{
			break;
		}

		// update statistics for adaptive candidate generation:
		if( true )
		{
			// calculate moments
			// update covariance matrix
		}

	}



}


/*
 * Calculate acceptance probability from Boltzmann distribution.
 * TODO: Should cost be capped at 1?
 */
/*
real
simulatedAnnealingModule::calculateAcceptanceProbability()
{
	real cap = 1.0;
	return std::min(std::exp( (candEnergy_ - crntEnergy_)/temp_ ), cap);
}
*/


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
 * Generates a candidate state in the neighbourhood of the current state.
 */
void simulatedAnnealingModule::generateCandidateState()
{
	// generate random direction in state space:
	real stateDir[stateDim_];
	for(int i = 0; i < stateDim_; i++)
	{
		stateDir[i] = candGenDistr_(rng_);
	}

	// copy current state to candidate state:
	cblas_scopy(stateDim_, crntState_, 1, candState_, 1);	

	// should adaptive candidate generation be used?
	if( useAdaptiveCandidateGeneration == true )
	{
		// candidate state is current state plus adapted random direction:
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, stateDim_, stateDim_, 1.0, 
//		            mat, 1, stateDir, 1, 1.0, candState_, 1);

	}
	else
	{
		// candidate state is current state plus random direction:
		cblas_saxpy(stateDim_, 1.0, stateDir, 1, candState_, 1); 
	}
}



/*
 * Generates a candidate state from the current state.
 */
/*
void
simulatedAnnealingModule::generateCandidate()
{

	int N = 3;
	double A[9] = {1, 2, 3, 2, 3, 4, 3, 4, 1};
	int LDA = 3;
	int IPIV[3];
	int INFO = 3;

	



	// loop over state dimension:
	for(int i=0; i<stateDim_; i++)
	{
		// new state is current state plus some small random offset:
		// TODO: introduce step length:
		candState_.at(i) = crntState_.at(i) + stepLengthFactor_*candGenDistr_(rng_);
	}
}
*/

