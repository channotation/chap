#include <algorithm>	// for std::max_element()
#include <cmath>		// for std::sqrt()

#include <gromacs/topology/atomprop.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>

#include "trajectoryAnalysis/trajectoryAnalysis.hpp"


using namespace gmx;



/*
 * Constructor for the trajectoryAnalysis class.
 */
trajectoryAnalysis::trajectoryAnalysis()
    : cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}



/*
 *
 */
void
trajectoryAnalysis::initOptions(IOptionsContainer          *options,
                                TrajectoryAnalysisSettings *settings)
{
	// set help text:
	static const char *const desc[] = {
		"This is a first prototype for the CHAP tool.",
		"There is NO HELP, you are on your own!"
	};
    settings -> setHelpText(desc);

	// require the user to provide a topology file input:
    settings -> setFlag(TrajectoryAnalysisSettings::efRequireTop);

	// get (required) selection option for the reference group: 
	options -> addOption(SelectionOption("reference")
	                     .store(&refsel_).required()
		                 .description("Reference group that defines the channel (normally 'Protein'): "));

	// get (required) selection options for the small particle groups:
	options -> addOption(SelectionOption("select")
                         .storeVector(&sel_).required().multiValue()
	                     .description("Groups to calculate distances to"));

    // get (optional) selection option for the neighbourhood search cutoff:
    options -> addOption(DoubleOption("cutoff")
	                     .store(&cutoff_)
                         .description("Cutoff for distance calculation (0 = no cutoff)"));	
}




/*
 * 
 */
void
trajectoryAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                                 const TopologyInformation &top)
{
	// set cutoff distance for grid search as specified in user input:
	nb_.setCutoff(cutoff_);
	std::cout<<"Setting cutoff to: "<<cutoff_<<std::endl;

	// set number of columns in data set (one column per small particle type):
	data_.setColumnCount(0, sel_.size());

	
	
	// load full topology:
	t_topology *topol = top.topology();	

	// access list of all atoms:
	t_atoms atoms = topol -> atoms;

	// create vector of van der Waals radii and allocate memory:
	vdwRadii.reserve(atoms.nr);

	// create atomprop struct:
	gmx_atomprop_t aps = gmx_atomprop_init();

	// loop over all atoms in system and get vdW-radii:
	for(int i=0; i<atoms.nr; i++)
	{
		real vdwRadius;

		// query vdW radius of current atom:
		if(gmx_atomprop_query(aps, 
		                      epropVDW, 
							  *(atoms.resinfo[atoms.atom[i].resind].name),
							  *(atoms.atomname[i]), &vdwRadius)) 
		{
			// TODO: include scale factor here?
		}
		else
		{
			// could not find vdW radius
			// TODO: handle this case
		}

		// add radius to vector of radii:
		vdwRadii.push_back(vdwRadius);
	}

	// delete atomprop struct:
	gmx_atomprop_destroy(aps);

	// find largest van der Waals radius in system:
	maxVdwRadius = *std::max_element(vdwRadii.begin(), vdwRadii.end());
}




/*
 *
 */
void
trajectoryAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                 TrajectoryAnalysisModuleData *pdata)
{
	// get data handle for this frame:
	AnalysisDataHandle dh = pdata -> dataHandle(data_);

	// get thread-local selection of reference particles:
	const Selection &ref_selection = pdata -> parallelSelection(refsel_);

	// get data for frame number frnr into data handle:
    dh.startFrame(frnr, fr.time);



	// PORE FINDING AND RADIUS CALCULATION
	// ------------------------------------------------------------------------

	// initialise neighbourhood search:
	AnalysisNeighborhoodSearch nbSearch = nb_.initSearch(pbc, ref_selection);


	// parameters:
	real probeStep = 1.0;
	RVec channelVector(0.0, 0.0, 1.0);
	RVec initialProbePosition(0, 0, 1);
	int maxProbeIter = 1;


	// initialise probe position:
	std::vector<RVec> probePosition = {initialProbePosition};



	// loop for advancing probe position:
	int i = 0;
	while(i < maxProbeIter)
	{
	
		// TODO: implement radius maximum radius calculating function
		// TODO: implement SA for radius finding		
		
		std::cout<<std::endl;
		std::cout<<"----------------------------------------------------"<<std::endl;


		real voidRadius = maximiseVoidRadius(initialProbePosition,
		                                     channelVector,
											 pbc,
											 ref_selection);
		std::cout<<"--> accepted radius = "<<voidRadius<<std::endl;


		// update probe position along channel axis:
		probePosition.front()[0] = probePosition.front()[0] + probeStep*channelVector[0];
		probePosition.front()[1] = probePosition.front()[1] + probeStep*channelVector[1];
		probePosition.front()[2] = probePosition.front()[2] + probeStep*channelVector[2];

	
		// increment loop counter:
		i++;
	}




	// wrapup


	// ANALYSIS OF SMALL PARTICLE POSITIONS
	// ------------------------------------------------------------------------

/*
	// loop over small particle selections:
	for(size_t g = 0; g < sel_.size(); ++g)
	{
		// get thread-local selection of small particles:
		const Selection &sel = pdata -> parallelSelection(sel_[g]);

		// get number of particles in selection:
		int n_part = sel.posCount();
	

		for(int i = 0; i < n_part; i++)
		{
			
			SelectionPosition p = sel.position(i);
			ConstArrayRef<int> idx = p.atomIndices();


			std::cout<<"vdW radius = "<<vdwRadii.at(idx.front())<<std::endl;

		}

	}

*/





	// finish analysis of current frame:
    dh.finishFrame();
}



/*
 *
 */
void
trajectoryAnalysis::finishAnalysis(int /*nframes*/)
{

}




void
trajectoryAnalysis::writeOutput()
{
	
}


/*
 * Function to calculate the radius of a spherical void around a given center.
 */
real
trajectoryAnalysis::calculateVoidRadius(RVec centre, 
                                        t_pbc *pbc, 
										const Selection refSelection)
{
	// convert centre coordinate to analysis neighbourhood position:
	std::vector<RVec> centrePosition = {centre};	
	AnalysisNeighborhoodPositions centrePos(centrePosition);


	/* TODO: calculate cutoff distance dynamicallyi
	 * can test this by comparison to radius calculated without cutoff!
	// initialise neighbourhood search:
	AnalysisNeighborhoodSearch nbSearch = nb_.initSearch(pbc, refSelection);

	// find distance to closest reference atom:
	real minDist = nbSearch.minimumDistance(centrePos);

	// set cutoff distance:
	nb_.setCutoff(1.0);
	*/


	// initialise neighbourhood pair search:
	AnalysisNeighborhoodSearch nbSearch = nb_.initSearch(pbc, refSelection);
	AnalysisNeighborhoodPairSearch nbPairSearch = nbSearch.startPairSearch(centrePos);
	AnalysisNeighborhoodPair pair;

	real voidRadius = 54321; // TODO:  infiinity
	real pairDist;
	real referenceVdwRadius;

	// loop over all pairs:
	while( nbPairSearch.findNextPair(&pair) )
	{	
		// find distance between particles in pair:
		// TODO: square root can probably be removed/only drawn for final radius?
		pairDist = std::sqrt( pair.distance2() );

		// get van der Waals radius of reference atom:
		referenceVdwRadius = vdwRadii.at(pair.refIndex());


		// update void radius if newly found distance is smaller:
		if( voidRadius > (pairDist - referenceVdwRadius) )
		{
			voidRadius = pairDist - referenceVdwRadius;
			std::cout<<"  pairDist = "<<pairDist<<"  vdW = "<<referenceVdwRadius<<"  smaller radius = "<<voidRadius<<std::endl;
		}
	}

	// return void radius:
	return voidRadius;
}


/*
 * Maximise the radius of a spherical void by relocation of void centre:
 */
real
trajectoryAnalysis::maximiseVoidRadius(RVec &centre,
									   RVec chanVec,
									   t_pbc *pbc,
									   const Selection refSelection)
{
	// parameters:
	int maxSimAnIter = 10;
	real initialTemperature = 100;
	real tempReductionFactor = 0.9;


	//
	RVec candidateCentre;
	real candidateVoidRadius;
	real voidRadius = -999999;
	real temp = initialTemperature;

	// generate random 3-vector:
	int seed=15011991;
	DefaultRandomEngine rng(seed);
    UniformRealDistribution<real> candGenDistr(-sqrt(3), sqrt(3)); // TODO: replace square roots with value for efficiency
	UniformRealDistribution<real> candAccDistr(0.0, 1.0);


	// simulated annealing iteration:
	for(int i=0; i<maxSimAnIter; i++)
	{
		// generate random 3-vector:
		RVec randVec(candGenDistr(rng), candGenDistr(rng), candGenDistr(rng));

		// remove components in the direction of channel vector:
		real scalarProduct = std::sqrt( randVec[0]*chanVec[0] + randVec[1]*chanVec[1] + randVec[2]*chanVec[2] );
		randVec[0] = randVec[0] - scalarProduct*chanVec[0];
		randVec[1] = randVec[1] - scalarProduct*chanVec[1];
		randVec[2] = randVec[2] - scalarProduct*chanVec[2];

		// calculate candidate for new centre position:
		candidateCentre = centre;
		candidateCentre[0] += randVec[0];
		candidateCentre[1] += randVec[1];
		candidateCentre[2] += randVec[2];

		// calculate radius of void around new centre:
		candidateVoidRadius = calculateVoidRadius(candidateCentre,
												  pbc,
												  refSelection);

		// calculate acceptance probability:
		// TODO: do we have to limit this to 1.0?
		real accProb = std::exp( (candidateVoidRadius - voidRadius)/temp );

		// accept move?
		if( accProb > candAccDistr(rng) )
		{
			// update centre position and void radius:
			centre = candidateCentre;
			voidRadius = candidateVoidRadius;
		}

		// reduce temperature:
		std::cout<<"temp = "<<temp<<"  x = "<<centre[0]<<"  y = "<<centre[1]<<"  z = "<<centre[2]<<"  radius = "<<voidRadius<<std::endl;
		temp *= tempReductionFactor;
	}


	// return maximised void radius:
	return voidRadius;
}






