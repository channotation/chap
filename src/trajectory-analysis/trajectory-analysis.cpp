#include <algorithm>	// for std::max_element()
#include <cmath>		// for std::sqrt()
#include <fstream>
#include <iomanip>
#include <string>
#include <ctime>

#include <gromacs/topology/atomprop.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>
#include <gromacs/fileio/confio.h>

#include "trajectory-analysis/trajectory-analysis.hpp"

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_1D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"

#include "io/molecular_path_obj_exporter.hpp"

#include "trajectory-analysis/simulated_annealing_module.hpp"
#include "trajectory-analysis/path_finding_module.hpp"
#include "trajectory-analysis/analysis_data_long_format_plot_module.hpp"
#include "trajectory-analysis/analysis_data_pdb_plot_module.hpp"

#include "path-finding/inplane_optimised_probe_path_finder.hpp"
#include "path-finding/optimised_direction_probe_path_finder.hpp"
#include "path-finding/naive_cylindrical_path_finder.hpp"

using namespace gmx;



/*
 * Constructor for the trajectoryAnalysis class.
 */
trajectoryAnalysis::trajectoryAnalysis()
    : cutoff_(0.0)
    , pfMethod_("inplane-optim")
    , pfProbeStepLength_(0.1)
    , pfProbeRadius_(0.0)
    , pfMaxFreeDist_(1.0)
    , pfMaxProbeSteps_(1e3)
    , pfInitProbePos_(3)
    , pfChanDirVec_(3)
    , saRandomSeed_(15011991)
    , saMaxCoolingIter_(1e3)
    , saNumCostSamples_(50)
    , saConvRelTol_(1e-10)
    , saInitTemp_(10.0)
    , saCoolingFactor_(0.99)
    , saStepLengthFactor_(0.01)
    , saUseAdaptiveCandGen_(false)
{
    //
    registerAnalysisDataset(&data_, "somedata");
    data_.setMultipoint(true);              // mutliple support points 
    
       // register dataset:
    registerAnalysisDataset(&dataResMapping_, "resMapping");
 



    // default initial probe position and chanell direction:
    pfInitProbePos_ = {0.0, 0.0, 0.0};
    pfChanDirVec_ = {0.0, 0.0, 1.0};


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

    // hardcoded defaults for multivalue options:
    std::vector<real> chanDirVec_ = {0.0, 0.0, 1.0};


	// require the user to provide a topology file input:
    settings -> setFlag(TrajectoryAnalysisSettings::efRequireTop);

	// get (required) selection option for the reference group: 
	options -> addOption(SelectionOption("reference")
	                     .store(&refsel_).required()
		                 .description("Reference group that defines the channel (normally 'Protein'): "));

	// get (required) selection options for the small particle groups:
	options -> addOption(SelectionOption("select")
                         .storeVector(&sel_).required()
	                     .description("Group of small particles to calculate density of (normally 'Water'):"));

   	// get (optional) selection options for the initial probe position selection:
	options -> addOption(SelectionOption("ippsel")
                         .store(&ippsel_)
                         .storeIsSet(&ippselIsSet_)
	                     .description("Reference group from which to determine the initial probe position for the pore finding algorithm. If unspecified, this defaults to the overall pore forming group. Will be overridden if init-probe-pos is set explicitly."));
 

    // get (optional) selection option for the neighbourhood search cutoff:
    options -> addOption(DoubleOption("cutoff")
	                     .store(&cutoff_)
                         .description("Cutoff for distance calculation (0 = no cutoff)"));


    // output options:
    options -> addOption(StringOption("ppfn")
                         .store(&poreParticleFileName_)
                         .defaultValue("pore_particles.dat")
                         .description("Name of file containing pore particle positions over time."));
    options -> addOption(StringOption("spfn")
                         .store(&smallParticleFileName_)
                         .defaultValue("small_particles.dat")
                         .description("Name of file containing small particle positions (i.e. water particle positions) over time."));
    options -> addOption(StringOption("o")
                         .store(&poreProfileFileName_)
                         .defaultValue("pore_profile.dat")
                         .description("Name of file containing pore radius, small particle density, and small particle energy as function of the permeation coordinate."));
    options -> addOption(IntegerOption("num-out-pts")
                         .store(&nOutPoints_)
                         .defaultValue(1000)
                         .description("Number of sample points of pore centre line that are written to output."));



    // get parameters of path-finding agorithm:
    options -> addOption(StringOption("pf-method")
                         .store(&pfMethod_)
                         .defaultValue("inplane-optim")
                         .description("Path finding method. Only inplane-optim is implemented so far."));
    options -> addOption(RealOption("probe-step")
                         .store(&pfProbeStepLength_)
                         .defaultValue(0.025)
                         .description("Step length for probe movement. Defaults to 0.025 nm."));
    options -> addOption(RealOption("probe-radius")
                         .store(&pfProbeRadius_)
                         .defaultValue(0.0)
                         .description("Radius of probe. Defaults to 0.0, buggy for other values!"));
    options -> addOption(RealOption("max-free-dist")
                         .store(&pfMaxFreeDist_)
                         .defaultValue(1.0)
                         .description("Maximum radius of pore. Defaults to 1.0, buggy for larger values."));
    options -> addOption(IntegerOption("max-probe-steps")
                         .store(&pfMaxProbeSteps_)
                         .description("Maximum number of steps the probe is moved in either direction."));
    options -> addOption(RealOption("init-probe-pos")
                         .storeVector(&pfInitProbePos_)
                         .storeIsSet(&pfInitProbePosIsSet_)
                         .valueCount(3)
                         .description("Initial position of probe in probe-based pore finding algorithms. If this is set explicitly, it will overwrite the COM-based initial position set with the ippselflag."));
    options -> addOption(RealOption("chan-dir-vec")
                         .storeVector(&pfChanDirVec_)
                         .storeIsSet(&pfChanDirVecIsSet_)
                         .valueCount(3)
                         .description("Channel direction vector; will be normalised to unit vector internally. Defaults to (0, 0, 1)."));
    options -> addOption(IntegerOption("sa-random-seed")
                         .store(&saRandomSeed_)
                         .required()
                         .description("Seed for RNG used in simulated annealing."));
    options -> addOption(IntegerOption("sa-max-cool")
                          .store(&saMaxCoolingIter_)
                          .defaultValue(1000)
                          .description("Maximum number of cooling iterations in one simulated annealing run. Defaults to 1000."));
    options -> addOption(IntegerOption("sa-cost-samples")
                         .store(&saNumCostSamples_)
                         .defaultValue(10)
                         .description("NOT IMPLEMENTED! Number of cost samples considered for convergence tolerance. Defaults to 10."));
    options -> addOption(RealOption("sa-conv-tol")
                         .store(&saConvRelTol_)
                         .defaultValue(1e-3)
                         .description("Relative tolerance for simulated annealing."));
    options -> addOption(RealOption("sa-init-temp")
                         .store(&saInitTemp_)
                         .defaultValue(0.1)
                         .description("Initital temperature for simulated annealing. Defaults to 0.1."));
    options -> addOption(RealOption("sa-cooling-fac")
                         .store(&saCoolingFactor_)
                         .defaultValue(0.98)
                         .description("Cooling factor using in simulated annealing. Defaults to 0.98."));
    options -> addOption(RealOption("sa-step")
                         .store(&saStepLengthFactor_)
                         .defaultValue(0.001)
                         .description("Step length factor used in candidate generation. Defaults to 0.001.")) ;
    options -> addOption(BooleanOption("debug-output")
                         .store(&debug_output_)
                         .description("When this flag is set, the program will write additional information.")) ;
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

	// prepare data container:
	data_.setDataSetCount(1);               // one data set for water   
    data_.setColumnCount(0, 5);             // x y z s r

    // add plot module to analysis data:
    int i = 2;
    AnalysisDataLongFormatPlotModulePointer lfplotm(new AnalysisDataLongFormatPlotModule(i));
    const char *poreParticleFileName = poreParticleFileName_.c_str();
    lfplotm -> setFileName(poreParticleFileName);
    lfplotm -> setPrecision(3);
    std::vector<char*> header = {"t", "x", "y", "z", "s", "r"};
    lfplotm -> setHeader(header);
    data_.addModule(lfplotm);  

    // add pdb plot module to analysis data:
    AnalysisDataPdbPlotModulePointer pdbplotm(new AnalysisDataPdbPlotModule(i));
    pdbplotm -> setFileName(poreParticleFileName);
    data_.addModule(pdbplotm);


    // RESIDUE MAPPING DATA
    //-------------------------------------------------------------------------


    // set dataset properties:
    dataResMapping_.setDataSetCount(1);
    dataResMapping_.setColumnCount(0, 4);   // refID s rho phi 
    dataResMapping_.setMultipoint(true);

    // add long format plot module:
    int j = 1;
    AnalysisDataLongFormatPlotModulePointer lfpltResMapping(new AnalysisDataLongFormatPlotModule(j));
    const char *fnResMapping = "res_mapping.dat";
    std::vector<char*> headerResMapping = {"t", "refId", "s", "rho", "phi"};
    lfpltResMapping -> setFileName(fnResMapping);
    lfpltResMapping -> setHeader(headerResMapping);
    lfpltResMapping -> setPrecision(15);    // TODO: different treatment for integers?
    dataResMapping_.addModule(lfpltResMapping);



    // PREPARE SELECTIONS FOR MAPPING
    //-------------------------------------------------------------------------

    gmx::SelectionCollection poreComCollection; 
    poreComCollection.setReferencePosType("res_com");
    poreComCollection.setOutputPosType("res_com");

    gmx::Selection test = poreComCollection.parseFromString("resname SOL")[0];

    poreComCollection.setTopology(top.topology(), 0);
    poreComCollection.compile();


    std::cout<<std::endl<<std::endl;
    std::cout<<"atomCount = "<<test.atomCount()<<"  "
             <<"posCount = "<<test.posCount()<<"  "
             <<std::endl;
    std::cout<<std::endl<<std::endl;

    


    
    // GET ATOM RADII FROM TOPOLOGY
    //-------------------------------------------------------------------------
	
	// load full topology:
	t_topology *topol = top.topology();	

	// access list of all atoms:
	t_atoms atoms = topol -> atoms;

	// create vector of van der Waals radii and allocate memory:
	vdwRadii_.reserve(atoms.nr);

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
		vdwRadii_.push_back(vdwRadius);
	}

	// delete atomprop struct:
	gmx_atomprop_destroy(aps);

	// find largest van der Waals radius in system:
	maxVdwRadius_ = *std::max_element(vdwRadii_.begin(), vdwRadii_.end());
}




/*
 *
 */
void
trajectoryAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                 TrajectoryAnalysisModuleData *pdata)
{
	// get data handles for this frame:
	AnalysisDataHandle dh = pdata -> dataHandle(data_);
    AnalysisDataHandle dhResMapping = pdata -> dataHandle(dataResMapping_);

	// get thread-local selection of reference particles:
	const Selection &refSelection = pdata -> parallelSelection(refsel_);

	// get data for frame number frnr into data handle:
    dh.startFrame(frnr, fr.time);
    dhResMapping.startFrame(frnr, fr.time);


    // UPDATE INITIAL PROBE POSITION FOR THIS FRAME
    //-------------------------------------------------------------------------

    // recalculate initial probe position based on reference group COM:
    if( pfInitProbePosIsSet_ == false )
    {  
        // helper variable for conditional assignment of selection:
        Selection tmpsel;
  
        // has a group for specifying initial probe position been set?
        if( ippselIsSet_ == true )
        {
            // use explicitly given selection:
            tmpsel = ippsel_;
        }
        else 
        {
            // default to overall group of pore forming particles:
            tmpsel = refsel_;
        }
     
        // load data into initial position selection:
        const gmx::Selection &initPosSelection = pdata -> parallelSelection(tmpsel);
 
        // initialse total mass and COM vector:
        real totalMass = 0.0;
        gmx::RVec centreOfMass(0.0, 0.0, 0.0);
        
        // loop over all atoms: 
        for(int i = 0; i < initPosSelection.atomCount(); i++)
        {
            // get i-th atom position:
            gmx::SelectionPosition atom = initPosSelection.position(i);

            // add to total mass:
            totalMass += atom.mass();

            // add to COM vector:
            // TODO: implement separate centre of geometry and centre of mass 
            centreOfMass[0] += atom.mass() * atom.x()[0];
            centreOfMass[1] += atom.mass() * atom.x()[1];
            centreOfMass[2] += atom.mass() * atom.x()[2];
        }

        // scale COM vector by total MASS:
        centreOfMass[0] /= 1.0 * totalMass;
        centreOfMass[1] /= 1.0 * totalMass;
        centreOfMass[2] /= 1.0 * totalMass; 

        // set initial probe position:
        pfInitProbePos_[0] = centreOfMass[0];
        pfInitProbePos_[1] = centreOfMass[1];
        pfInitProbePos_[2] = centreOfMass[2];
    }

    // inform user:
    if( debug_output_ == true )
    {
        std::cout<<std::endl
                 <<"Initial probe position for this frame is:  "
                 <<pfInitProbePos_[0]<<", "
                 <<pfInitProbePos_[1]<<", "
                 <<pfInitProbePos_[2]<<". "
                 <<std::endl;
    }


    // GET VDW RADII FOR SELECTION
    //-------------------------------------------------------------------------
    // TODO: Move this to separate class and test!
    // TODO: Should then also work for coarse-grained situations!

	// create vector of van der Waals radii and allocate memory:
    std::vector<real> selVdwRadii;
	selVdwRadii.reserve(refSelection.atomCount());
    std::cout<<"selVdwRadii.size() = "<<selVdwRadii.size()<<std::endl;

	// loop over all atoms in system and get vdW-radii:
	for(int i=0; i<refSelection.atomCount(); i++)
    {
        // get global index of i-th atom in selection:
        gmx::SelectionPosition atom = refSelection.position(i);
        int idx = atom.mappedId();

		// add radius to vector of radii:
		selVdwRadii.push_back(vdwRadii_[idx]);
	}


	// PORE FINDING AND RADIUS CALCULATION
	// ------------------------------------------------------------------------

	// initialise neighbourhood search:
	AnalysisNeighborhoodSearch nbSearch = nb_.initSearch(pbc, refSelection);

   /* 
    std::cout<<"pfMethod = "<<pfMethod_<<std::endl
             <<"pfProbeStepLength = "<<pfProbeStepLength_<<std::endl
             <<"pfProbeRadius = "<<pfProbeRadius_<<std::endl
             <<"pfMaxFreeDist = "<<pfMaxFreeDist_<<std::endl
             <<"pfMaxProbeSteps = "<<pfMaxProbeSteps_<<std::endl
             <<"pfInitProbePos = "<<pfInitProbePos_[0]<<"  "
                                  <<pfInitProbePos_[1]<<"  "
                                  <<pfInitProbePos_[2]<<std::endl
             <<"pfChanDirVec = "<<pfChanDirVec_[0]<<"  "
                                <<pfChanDirVec_[1]<<"  "
                                <<pfChanDirVec_[2]<<std::endl
             <<"saRandomSeed = "<<saRandomSeed_<<std::endl
             <<"saMaxCoolingIter = "<<saMaxCoolingIter_<<std::endl
             <<"saNumCostSamples = "<<saNumCostSamples_<<std::endl
             <<"saXi = "<<saXi_<<std::endl
             <<"saConvRelTol = "<<saConvRelTol_<<std::endl
             <<"saInitTemp = "<<saInitTemp_<<std::endl
             <<"saCoolingFactor = "<<saCoolingFactor_<<std::endl
             <<"saStepLengthFactor = "<<saStepLengthFactor_<<std::endl
             <<"saUseAdaptiveCandGen = "<<saUseAdaptiveCandGen_<<std::endl;
    */


    // create path finding module:
    std::unique_ptr<AbstractPathFinder> pfm;
    if( pfMethod_ == "inplane-optim" )
    {
    	RVec initProbePos(pfInitProbePos_[0], pfInitProbePos_[1], pfInitProbePos_[2]);
    	RVec chanDirVec(pfChanDirVec_[0], pfChanDirVec_[1], pfChanDirVec_[2]);
        pfm.reset(new InplaneOptimisedProbePathFinder(pfProbeStepLength_, pfProbeRadius_, 
                                            pfMaxFreeDist_, pfMaxProbeSteps_, 
                                            initProbePos, chanDirVec, selVdwRadii, 
                                            &nbSearch, saRandomSeed_, 
                                            saMaxCoolingIter_, saNumCostSamples_, 
                                            saXi_, saConvRelTol_, saInitTemp_, 
                                            saCoolingFactor_, saStepLengthFactor_, 
                                            saUseAdaptiveCandGen_));
    }
    else if( pfMethod_ == "optim-direction" )
    {
        std::cout<<"OPTIM-DIRECTION"<<std::endl;

    	RVec initProbePos(pfInitProbePos_[0], pfInitProbePos_[1], pfInitProbePos_[2]);
    	RVec chanDirVec(pfChanDirVec_[0], pfChanDirVec_[1], pfChanDirVec_[2]);
        pfm.reset(new OptimisedDirectionProbePathFinder(pfProbeStepLength_, pfProbeRadius_, 
                                            pfMaxFreeDist_, pfMaxProbeSteps_, 
                                            initProbePos, selVdwRadii, 
                                            &nbSearch, saRandomSeed_, 
                                            saMaxCoolingIter_, saNumCostSamples_, 
                                            saXi_, saConvRelTol_, saInitTemp_, 
                                            saCoolingFactor_, saStepLengthFactor_, 
                                            saUseAdaptiveCandGen_));
    }   
    else if( pfMethod_ == "naive-cylindrical" )
    {
    	RVec initProbePos(pfInitProbePos_[0], pfInitProbePos_[1], pfInitProbePos_[2]);
    	RVec chanDirVec(pfChanDirVec_[0], pfChanDirVec_[1], pfChanDirVec_[2]);
        pfm.reset(new NaiveCylindricalPathFinder(pfProbeStepLength_,
                                                 pfMaxProbeSteps_,
                                                 pfMaxFreeDist_,
                                                 initProbePos,
                                                 chanDirVec));        
    }










    const gmx::Selection &test = pdata -> parallelSelection(refsel_);


//    std::cout<<"atomCount = "<<test.atomCount()<<"  "
//             <<std::endl;




























//    std::cout<<"================================================="<<std::endl;
//    std::cout<<std::endl;
//    std::cout<<std::endl;


    // run path finding algorithm on current frame:
//    std::cout<<"finding permeation pathway ... ";
    clock_t tPathFinding = std::clock();
    pfm -> findPath();
    tPathFinding = (std::clock() - tPathFinding)/CLOCKS_PER_SEC;
//    std::cout<<"done in  "<<tPathFinding<<" sec"<<std::endl;


    // retrieve molecular path object:
//    std::cout<<"preparing pathway object ... ";
    clock_t tMolPath = std::clock();
    MolecularPath molPath = pfm -> getMolecularPath();
    tMolPath = (std::clock() - tMolPath)/CLOCKS_PER_SEC;
//    std::cout<<"done in  "<<tMolPath<<" sec"<<std::endl;

    
    // map residues onto pathway:
//    std::cout<<"mapping residues onto pathway ... ";
    clock_t tMapRes = std::clock();
//    const gmx::Selection &refResidueSelection = pdata -> parallelSelection(refsel_);
//    std::map<int, gmx::RVec> mappedCoords = molPath.mapSelection(refResidueSelection, pbc);
    tMapRes = (std::clock() - tMapRes)/CLOCKS_PER_SEC;
//    std::cout<<"done in  "<<tMapRes<<" sec"<<std::endl;

 //   std::cout<<mappedCoords.size()<<" particles have been mapped"<<std::endl;


    // check if points lie inside pore:
//    std::cout<<"checking if particles are inside pore ... ";
//    clock_t tCheckInside = std::clock();
//    std::map<int, bool> isInside = molPath.checkIfInside(mappedCoords);
//    tCheckInside = (std::clock() - tCheckInside)/CLOCKS_PER_SEC;
//    std::cout<<"done in  "<<tCheckInside<<" sec"<<std::endl;

/*
    for(std::map<int, gmx::RVec>::iterator it = mappedCoords.begin(); it != mappedCoords.end(); it++)
    {
         dhResMapping.setPoint(0, it -> first);         // refId
         dhResMapping.setPoint(1, it -> second[0]);     // s
         dhResMapping.setPoint(2, it -> second[1]);     // rho
         dhResMapping.setPoint(3, it -> second[3]);     // phi
         dhResMapping.finishPointSet();
    }
*/
    

    

    std::cout<<std::endl;
    std::cout<<std::endl;
//    std::cout<<"================================================="<<std::endl;


    // reset smart pointer:
    // TODO: Is this necessary and parallel compatible?
    pfm.reset();





    // ADD PATH DATA TO PARALLELISABLE CONTAINER
    //-------------------------------------------------------------------------
/*
    // access path finding module result:
    real extrapDist = 1.0;
    std::vector<real> arcLengthSample = molPath.sampleArcLength(nOutPoints_, extrapDist);
    std::vector<gmx::RVec> pointSample = molPath.samplePoints(arcLengthSample);
    std::vector<real> radiusSample = molPath.sampleRadii(arcLengthSample);
 
    // loop over all support points of path:
    for(int i = 0; i < nOutPoints_; i++)
    {
        // add to container:
        dh.setPoint(0, pointSample[i][0]);     // x
        dh.setPoint(1, pointSample[i][1]);     // y
        dh.setPoint(2, pointSample[i][2]);     // z
        dh.setPoint(3, arcLengthSample[i]);    // s
        dh.setPoint(4, radiusSample[i]);       // r

        dh.finishPointSet(); 
    }
  */

    // WRITE PORE TO OBJ FILE
    //-------------------------------------------------------------------------


    MolecularPathObjExporter molPathExp;
    molPathExp("pore.obj",
               molPath);


    // FINISH FRAME
    //-------------------------------------------------------------------------

	// finish analysis of current frame:
    dh.finishFrame();
    dhResMapping.finishFrame();
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
    std::cout<<"datSetCount = "<<data_.dataSetCount()<<std::endl
             <<"columnCount = "<<data_.columnCount()<<std::endl
             <<"frameCount = "<<data_.frameCount()<<std::endl
             <<std::endl;
}


