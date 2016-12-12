#include <algorithm>	// for std::max_element()
#include <cmath>		// for std::sqrt()
#include <fstream>
#include <iomanip>
#include <string>

#include <gromacs/topology/atomprop.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>

#include "trajectoryAnalysis/trajectoryAnalysis.hpp"

#include "trajectoryAnalysis/simulated_annealing_module.hpp"
#include "trajectoryAnalysis/path_finding_module.hpp"

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
	const Selection &refSelection = pdata -> parallelSelection(refsel_);


    std::cout<<"atomCount = "<<refSelection.atomCount()<<std::endl;



	// get data for frame number frnr into data handle:
    dh.startFrame(frnr, fr.time);



	// PORE FINDING AND RADIUS CALCULATION
	// ------------------------------------------------------------------------

	// initialise neighbourhood search:
	AnalysisNeighborhoodSearch nbSearch = nb_.initSearch(pbc, refSelection);


	// parameters:
	real probeStep = 0.5;
	RVec channelVector(0.0, 0.0, 1.0);

    
    
	RVec initialProbePosition(55, 55, 40.0);
    int maxProbeIter = 1;

    std::cout<<"blah"<<std::endl;
    PathFindingModule pfm(initialProbePosition, channelVector, nbSearch, vdwRadii);
    pfm.findPath();


    std::vector<gmx::RVec> path = pfm.getPath();
    std::vector<real> radii = pfm.getRadii();


    // write path to PDB file:
    std::fstream pdbfile;
    std::string pdbfilename = "test.pdb"; 
    pdbfile.open(pdbfilename.c_str(), std::fstream::out);
    pdbfile.precision(3);

    pdbfile<<"HEADER"<<"   test"<<std::endl;
    pdbfile<<"TITLE"<<"   test"<<std::endl;
   

    for(int i=0; i<path.size(); i++)
    {
        pdbfile<<std::setw(6)<<"ATOM  "                 // record name
            <<std::setw(5)<<i+1                    // atom serial number (one-based)
            <<std::setw(1)<<" "
            <<std::setw(4)<<"PORE"                 // atom name
            <<std::setw(1)<<" "                    // alternate location indicator
            <<std::setw(3)<<"POR"                  // residue name
            <<std::setw(1)<<""
            <<std::setw(1)<<"X"                    // chain identifier
            <<std::setw(4)<<"000"                  // residue sequence number
            <<std::setw(1)<<" "                    // code for insertion of residues
            <<std::setw(3)<<""
            <<std::setw(8)<<path[i][0]             // x
            <<std::setw(8)<<path[i][1]             // y
            <<std::setw(8)<<path[i][2]             // z
            <<std::setw(6)<<radii[i]                  // occupancy
            <<std::setw(6)<<radii[i]               // temperature factor
            <<std::setw(10)<<"" 
            <<std::setw(2)<<"XX"                   // element symbol
            <<std::setw(2)<<0                      // charge
            <<std::endl;

    }

    pdbfile<<"  10.34380  10.34380  10.83500"<<std::endl;

    pdbfile.close();



    // write path to GRO file:
    std::fstream file;
    std::string filename = "test.gro"; 
    file.open(filename.c_str(), std::fstream::out);
    file.precision(3);

    file<<"titlestring"<<std::endl;
    file<<radii.size()<<std::endl;

    for(int i=0; i<path.size(); i++)
    {
        file<<std::setw(5)<<1                 // residue number
            <<std::setw(5)<<"PORE"            // residue name
            <<std::setw(5)<<"PORE"            // atom name
            <<std::setw(5)<<i                 // atom number
            <<std::setw(8)<<path[i][0]        // x
            <<std::setw(8)<<path[i][1]        // y
            <<std::setw(8)<<path[i][2]        // z
            <<std::setw(8)<<0.000             // vx
            <<std::setw(8)<<0.000             // vy
            <<std::setw(8)<<0.000             // vz
            <<std::endl;

    }

    file<<"  10.34380  10.34380  10.83500"<<std::endl;

    file.close();



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
//			std::cout<<"  pairDist = "<<pairDist<<"  vdW = "<<referenceVdwRadius<<"  smaller radius = "<<voidRadius<<std::endl;
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
	int maxSimAnIter = 1000;
	real initialTemperature = 100;
	real tempReductionFactor = 0.99;


	//
	RVec candidateCentre;
	real candidateVoidRadius;
	real voidRadius = -999999;
	real bestVoidRadius = voidRadius;
	real temp = initialTemperature;
	

	// generate random 3-vector:
	int seed=15011991;
	real candStepLength = 0.01;
	DefaultRandomEngine rng(seed);
    UniformRealDistribution<real> candGenDistr(-candStepLength*sqrt(3), candStepLength*sqrt(3)); // TODO: replace square roots with value for efficiency
	UniformRealDistribution<real> candAccDistr(0.0, 1.0);


	// simulated annealing iteration:
	for(int i=0; i<maxSimAnIter; i++)
	{
		// generate random 3-vector:
		RVec randVec(candGenDistr(rng), candGenDistr(rng), candGenDistr(rng));

		// remove components in the direction of channel vector:
		real scalarProduct = randVec[0]*chanVec[0] + randVec[1]*chanVec[1] + randVec[2]*chanVec[2];
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

//			std::cout<<"accepted!"<<std::endl;
//
			if( voidRadius > bestVoidRadius )
			{
				bestVoidRadius = voidRadius;
			}
		}


		// reduce temperature:
		std::cout<<"temp = "<<temp<<"  x = "<<centre[0]<<"  y = "<<centre[1]<<"  z = "<<centre[2]<<"  radius = "<<voidRadius<<std::endl;
		temp *= tempReductionFactor;
	}


	// return maximised void radius:
	return bestVoidRadius;
}






