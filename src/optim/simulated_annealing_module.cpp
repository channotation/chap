#include <iostream>
#include <numeric>
#include <functional>

#include <cblas.h>
#include <lapacke.h>

#include "optim/simulated_annealing_module.hpp"
#include "statistics/misc_stats_utilities.hpp"


/*!
 * Simple constructor. Creates an SimulatedAnnealingModule object, but does
 * not set any of its properties.
 */
SimulatedAnnealingModule::SimulatedAnnealingModule()
{

}


/*!
 * Custom destructor to free memory.
 */
SimulatedAnnealingModule::~SimulatedAnnealingModule()
{
	// free memory occupied by internal state vectors:
	delete[] crntState_;
	delete[] candState_;
	delete[] bestState_;
}


/*!
 * Sets parameters of the simulated annealing algorithm. Will throw an error
 * if any required parameter without defaults is not set. Ignores unknown 
 * parameters.
 */
void
SimulatedAnnealingModule::setParams(std::map<std::string, real> params)
{
    // PRNG seed:
    if( params.find("saSeed") != params.end() )
    {
        seed_ = params["saSeed"];
    }
    else
    {
        // TODO: random seed!
    }
    
    // number of cooling iterations:
    if( params.find("saMaxCoolingIter") != params.end() )
    {
        maxCoolingIter_ = params["saMaxCoolingIter"];
    }
    else
    {
        std::cerr<<"ERROR: No maximum number of cooling iterations given!"<<std::endl;
        std::abort();
    }

    // initial temperature:
    if( params.find("saInitTemp") != params.end() )
    {
        temp_ = params["saInitTemp"];
    }
    else
    {
        std::cerr<<"ERROR: No initial temperature given!"<<std::endl;
        std::abort();
    }

    // cooling factor:
    if( params.find("saCoolingFactor") != params.end() )
    {
        coolingFactor_ = params["saCoolingFactor"];
    }
    else
    {
        std::cerr<<"ERROR: No cooling factor given!"<<std::endl;
        std::abort();
    }

    // step length factor:
    if( params.find("saStepLengthFactor") != params.end() )
    {
        stepLengthFactor_ = params["saStepLengthFactor"];
    }
    else
    {
        std::cerr<<"ERROR: No step length factor given!"<<std::endl;
        std::abort();
    }
}


/*!
 * Sets the objective function object.
 */
void
SimulatedAnnealingModule::setObjFun(ObjectiveFunction objFun)
{
     this -> objFun_ = objFun;
}


/*!
 * Sets the initial point in optimisation space from which simulated annealing
 * is started. This function also allocates the memory needed by the arrays
 * containing the internal state.
 */
void
SimulatedAnnealingModule::setInitGuess(std::vector<real> guess)
{
    // set optimisation space dimension:
    stateDim_ = guess.size();

    // allocate memory for internal state vectors:
	crntState_ = new real[stateDim_];
	candState_ = new real[stateDim_];
	bestState_ = new real[stateDim_];

	// initialise state vectors:
    std::copy(guess.begin(), guess.end(), crntState_);
    std::copy(guess.begin(), guess.end(), candState_);
    std::copy(guess.begin(), guess.end(), bestState_);
}


/*!
 * Implements the optimisation class interface. Wraps around the anneal() 
 * function, which stems from the original version of this class.
 */
void
SimulatedAnnealingModule::optimise()
{
    anneal();
}


/*!
 * Returns the optimisation result (i.e. the best point found) and the 
 * corresponding objective function value as am OptimSpacePoint.
 */
OptimSpacePoint
SimulatedAnnealingModule::getOptimPoint()
{
    OptimSpacePoint res;
    res.first.assign(bestState_, bestState_ + stateDim_);
    res.second = bestCost_;

    return res;
}


/*!
 * Public interface for the annealing function.
 *
 * This was intended to make the distinction between adaptive and isotropic 
 * annealing, but since the adaptive annealing procedure has been removed, all
 * this handles is an evaluation of the cost of the initial state before
 * calling annealIsotropic().
 */
void
SimulatedAnnealingModule::anneal()
{
   	// get cost of inital states:
    std::vector<real> crntStateVec(crntState_, crntState_ + stateDim_);
    std::vector<real> candStateVec(candState_, candState_ + stateDim_);
    std::vector<real> bestStateVec(bestState_, bestState_ + stateDim_);
    crntCost_ = objFun_(crntStateVec);
    candCost_ = objFun_(candStateVec);
    bestCost_ = objFun_(bestStateVec);

    // adaptive annealing not implemented:
    annealIsotropic();
}


/*!
 * Nonadaptive version of the annealing procedure. At each temperature, the 
 * cost function is evaluated exactly once and candidate states are always 
 * generated by making a small step in a isotropically random direction.
 */
void
SimulatedAnnealingModule::annealIsotropic()
{
	// initialise counter:
	int nCoolingIter = 0;

	// start annealing loop:
	while(true)
	{
        // generate a candidate state:
        generateCandidateStateIsotropic();

        // TODO: check boundary conditions!
        
        // evaluate cost function:
        std::vector<real> candStateVec;
        candStateVec.assign(candState_, candState_ + stateDim_);
        candCost_ = objFun_(candStateVec);

        // accept candidate?
        if( acceptCandidateState() == true )
        {
            // candidate state becomes current state:
            cblas_scopy(stateDim_, candState_, 1, crntState_, 1);
            crntCost_ = candCost_;
            // is new state also the best state?
            if( candCost_ > bestCost_ )
            {
                cblas_scopy(stateDim_, candState_, 1, bestState_, 1);
                bestCost_ = candCost_;
            }
        }

        // reduce temperature:
        cool();
        nCoolingIter++;

        // maximum step number reached?
        if( nCoolingIter >= maxCoolingIter_ )
        {
            return;
        }
	}
}


/*!
 * Reduces temperature of SA module. Currently only simple exponential 
 * cooling is implemented, i.e.
 *
 * \f[
 *      T_{i+1} = \gamma T_{i}
 * \f]
 *
 * where \f$ \gamma \in (0,1) \f$ is a cooling factor.
 */
void
SimulatedAnnealingModule::cool()
{
	temp_ *= coolingFactor_;
}


/*!
 * Generates a candidate state in the neighbourhood of the current state, where
 * the step direction generated isotropically at random.
 */
void 
SimulatedAnnealingModule::generateCandidateStateIsotropic()
{
	// generate random direction in state space:
	real stateDir[stateDim_];
	for(int i = 0; i < stateDim_; i++)
	{
		stateDir[i] = candGenDistr_(rng_);
	}

	// copy current state to candidate state:
	cblas_scopy(stateDim_, crntState_, 1, candState_, 1);	

	// candidate state is current state plus random direction:
	cblas_saxpy(stateDim_, stepLengthFactor_, stateDir, 1, candState_, 1);
}


/*!
 * Decides whether to accept or reject a candidate state. Returns true if 
 * candidate state is accepted.
 *
 * Acceptance probability is calculated as 
 *
 * \f[
 *      P(\text{accept}) = \min\left( \exp{ \frac{c_{\text{cand}} - 
 *      c_{\text{crnt}}}{T} }, 1 \right)
 * \f]
 *
 * where \f$ c_{\text{*}} \f$ is the candidate and current cost 
 * respectively and \f$ T \f$ is the current temperature. This is then
 * compared to a uniform random number on the interval \f$ [0,1) \f$ to
 * determine whether to accept a candidate state.
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

