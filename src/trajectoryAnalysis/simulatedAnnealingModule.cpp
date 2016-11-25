#include <iostream>

#include <cblas.h>
#include <lapacke.h>

#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"


/*
 * Constructor.
 *
 * Note that the candidate generation range is [-sqrt(3), sqrt(3)), which gives
 * a standard deviation of 1.0 for a uniform distribution.
 */
SimulatedAnnealingModule::SimulatedAnnealingModule(int stateDim,
                                                   int randomSeed,
												   int maxCoolingIter,
												   int numCovSamples,
												   real initTemp,
												   real coolingFactor,
												   real stepLengthFactor,
												   real *initState,
												   costFunction cf,
												   bool useAdaptiveCandidateGeneration)
	: stateDim_(stateDim)
	, seed_(randomSeed)
	, maxCoolingIter_(maxCoolingIter)
	, numCovSamples_(numCovSamples)
	, temp_(initTemp)
	, coolingFactor_(coolingFactor)
	, stepLengthFactor_(stepLengthFactor)
	, rng_(seed_)
	, candGenDistr_(-1.732051, 1.732051)
	, candAccDistr_(0.0, 1.0)
	, crntState_()
	, candState_()
	, bestState_()
	, evaluateCost_(cf)
	, useAdaptiveCandidateGeneration_(useAdaptiveCandidateGeneration)
{
	// allocate memory for internal state vectors:
	crntState_ = new real[stateDim_];
	candState_ = new real[stateDim_];
	bestState_ = new real[stateDim_];

	// initialise state vectors:
	cblas_scopy(stateDim_, initState, 1, crntState_, 1);
	cblas_scopy(stateDim_, initState, 1, candState_, 1);
	cblas_scopy(stateDim_, initState, 1, bestState_, 1);

	// get cost of inital states:
	crntCost_ = evaluateCost_(crntState_);
	candCost_ = evaluateCost_(candState_);
	bestCost_ = evaluateCost_(bestState_);

	// should adaptive candidate generation be used:
	if( useAdaptiveCandidateGeneration_ == true )
	{
		// allocate memory for covariance and adaptation matrix:
		stateSampleMatrix_ = new real[stateDim_ * numCovSamples_];
		covarianceMatrix_ = new real[stateDim_ * stateDim_];
		adaptationMatrix_ = new real[stateDim_ * stateDim_];
	}
}


/*
 * Custom destructor to free memory.
 */
SimulatedAnnealingModule::~SimulatedAnnealingModule()
{
	// free memory occupied by internal state vectors:
	delete[] crntState_;
	delete[] candState_;
	delete[] bestState_;


	if( useAdaptiveCandidateGeneration_ == true )
	{
		// free memory occupied by adaptation algorithm matrices:
		delete[] stateSampleMatrix_;
		delete[] covarianceMatrix_;
		delete[] adaptationMatrix_;
	}
}


/*
 * Perform annealing procedure.
 *
 * TODO: Handle constraints!
 */
void
SimulatedAnnealingModule::anneal()
{
	// initialise counter:
	nCoolingIter_ = 0;


	// start annealing loop:
	while(true)
	{
		// reset counter of accepted and rejected moves:
		int nAccepted = 0;
		int nRejected = 0;

		// 
		for(int i = 0; i < numCovSamples_; i++)
		{
			// generate a candidate state:
			generateCandidateState();

			// TODO: check boundary conditions!
			
			// evaluate cost function:
			candCost_ = evaluateCost_(candState_);


			// accept candidate?
			if( acceptCandidateState() == true )
			{
				// candidate state becomes current state:
				cblas_scopy(stateDim_, candState_, 1, crntState_, 1);
				crntCost_ = candCost_;

				// increment acceptance counter:
				nAccepted++;

				// is new state also the best state?
				if( candCost_ > bestCost_ )
				{
					cblas_scopy(stateDim_, candState_, 1, bestState_, 1);
					bestCost_ = candCost_;
				}

			}
			else
			{
				// increment rejection counter:
				nRejected++;
			}

			// copy current state to sample matrix:
			if( useAdaptiveCandidateGeneration_ )
			{
				LAPACKE_slacpy(LAPACK_ROW_MAJOR, 'A', stateDim_, 1, crntState_,
							   1, &stateSampleMatrix_[i], numCovSamples_);
			}
			
		}

		// reduce temperature:
		cool();
		nCoolingIter_++;
		
		// update statistics for adaptive candidate generation:
		if( useAdaptiveCandidateGeneration_ )
		{
			updateAdaptationMatrix();
		}
/*	
		std::cout<<"i = "<<nCoolingIter_
				 <<"  temp = "<<temp_
		         <<"  crntCost = "<<crntCost_
				 <<"  crntState = "<<crntState_[0]<<","<<crntState_[1]
				 <<"  candCost = "<<candCost_
				 <<"  candState = "<<candState_[0]<<","<<candState_[1]
				 <<"  bestCost = "<<bestCost_
				 <<"  bestState = "<<bestState_[0]<<","<<bestState_[1]
				 <<std::endl;
*/

		if( useAdaptiveCandidateGeneration_ )
		{

			for(int n=0; n<stateDim_; n++)
			{
				for(int m=0; m<numCovSamples_; m++)
				{
					std::cout<<stateSampleMatrix_[m + n*numCovSamples_]<<"  ";
				}
				std::cout<<std::endl;
			}
		}

		// maximum step number reached?
		if( nCoolingIter_ >= maxCoolingIter_ )
		{
			break;
		}

		// no candidates accepted?
		if( useAdaptiveCandidateGeneration_ && nAccepted <= 0 )
		{
			break;
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
SimulatedAnnealingModule::cool()
{
	temp_ *= coolingFactor_;
}



/*
 * Generates a candidate state in the neighbourhood of the current state.
 */
void 
SimulatedAnnealingModule::generateCandidateState()
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
//	if( useAdaptiveCandidateGeneration_ == true )
//	{
		// candidate state is current state plus adapted random direction:
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, stateDim_, stateDim_, 1.0, 
//		            adaptationMatrix_, 1, stateDir, 1, 1.0, candState_, 1);

//	}
//	else
//	{
		// candidate state is current state plus random direction:
		cblas_saxpy(stateDim_, stepLengthFactor_, stateDir, 1, candState_, 1); 
//	}
}


/*
 * Decides whether to accept or reject a candidate state. Returns true if 
 * candidate state is accepted.
 */
bool
SimulatedAnnealingModule::acceptCandidateState()
{
	// calculate acceptance probability according to Boltzmann statistics:
	real accProb = std::min(std::exp( (candCost_ - crntCost_)/temp_ ), 1.0f);

	// draw unfiform random number on interval [0,1):
	real r = candAccDistr_(rng_);
	
	// should candidate be accepted:
	return (r < accProb);
}


/*
 *
 */
void
SimulatedAnnealingModule::updateAdaptationMatrix()
{
	// calculate fist moment:
	real firstMoment[stateDim_] = {0};
	for(int i = 0; i < numCovSamples_; i++)
	{
		std::cout<<"i = "<<i<<"  m = "<<i<<"  n = "<<i + numCovSamples_<<"  state(i) = "<<stateSampleMatrix_[i]<<","<<stateSampleMatrix_[i + numCovSamples_];
		cblas_saxpy(stateDim_, 1.0f, &stateSampleMatrix_[i], numCovSamples_,
					firstMoment, 1);
		std::cout<<"  firstMoment = "<<firstMoment[0]<<","<<firstMoment[1]<<std::endl;
	}
	cblas_sscal(stateDim_, 1.0f/numCovSamples_, firstMoment,	1);
	std::cout<<"firstMoment = "<<firstMoment[0]<<","<<firstMoment[1]<<std::endl;



	// de-mean state sample matrix:
	for(int i = 0; i < numCovSamples_; i++)
	{
		cblas_saxpy(stateDim_, -1.0f, firstMoment, 1, &stateSampleMatrix_[i], numCovSamples_);
	}


	// calculate second moment:
	real secondMoment[stateDim_ * stateDim_] = {0};
	


	std::cout<<"---"<<std::endl;
	for(int i = 0; i < stateDim_; i++)
	{
		for(int j = 0; j < numCovSamples_; j++)
		{
			std::cout<<""<<stateSampleMatrix_[j + i*numCovSamples_]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"---"<<std::endl;



	for(int i = 0; i < numCovSamples_; i++)
	{
		cblas_sgemm(CblasRowMajor,
					CblasNoTrans,
		              CblasTrans,
					  stateDim_,
					  stateDim_,
					  numCovSamples_,
					  1.0f,
					  stateSampleMatrix_,
					  numCovSamples_,
					  stateSampleMatrix_,
					  numCovSamples_,
					  1.0f,
					  secondMoment,
					  stateDim_);		
	}

	
	cblas_sscal(stateDim_*stateDim_, 1.0f/numCovSamples_, secondMoment, 1);







	std::cout<<"---"<<std::endl;
	for(int i = 0; i < stateDim_; i++)
	{
		for(int j = 0; j < stateDim_; j++)
		{
			std::cout<<""<<stateSampleMatrix_[j + i*stateDim_]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"---"<<std::endl;
	
}


