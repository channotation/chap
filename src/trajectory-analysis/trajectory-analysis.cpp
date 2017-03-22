#include <algorithm>	// for std::max_element()
#include <cmath>		// for std::sqrt()
#include <fstream>
#include <iomanip>
#include <string>

#include <gromacs/topology/atomprop.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>
#include <gromacs/fileio/confio.h>

#include "trajectory-analysis/trajectory-analysis.hpp"

#include "trajectory-analysis/simulated_annealing_module.hpp"
#include "trajectory-analysis/path_finding_module.hpp"
#include "path-finding/inplane_optimised_probe_path_finder.hpp"
#include "path-finding/optimised_direction_probe_path_finder.hpp"


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
    registerAnalysisDataset(&data_, "avedist");


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

	// set number of columns in data set (one column per small particle type):
	data_.setColumnCount(0, sel_.size());

	
	
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


	// get data handle for this frame:
	AnalysisDataHandle dh = pdata -> dataHandle(data_);

	// get thread-local selection of reference particles:
	const Selection &refSelection = pdata -> parallelSelection(refsel_);

	// get data for frame number frnr into data handle:
    dh.startFrame(frnr, fr.time);


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
        real max_z = 0;
        real min_z = 0;
        // loop over all atoms: 
        for(int i = 0; i < initPosSelection.atomCount(); i++)
        {
            // get i-th atom position:
            gmx::SelectionPosition atom = initPosSelection.position(i);

            std::cout<<"nAtoms = "<<atom.atomCount()<<", "
                     <<"x = "<<atom.x()[0]<<", "
                     <<"y = "<<atom.x()[1]<<", "
                     <<"z = "<<atom.x()[2]<<", "
                     <<"mass = "<<atom.mass()
                     <<std::endl;
 
            // add to total mass:
            totalMass += atom.mass();

            // add to COM vector:
//`            centreOfMass[0] += atom.mass() * atom.x()[0];
//            centreOfMass[1] += atom.mass() * atom.x()[1];
//            centreOfMass[2] += atom.mass() * atom.x()[2];
            centreOfMass[0] += atom.x()[0];
            centreOfMass[1] += atom.x()[1];
            centreOfMass[2] += atom.x()[2];

            if( atom.x()[2] > max_z )
            {
                max_z = atom.x()[2];
            }
            if( atom.x()[2] < min_z )
            {
                min_z = atom.x()[2];
            }

        }

        std::cout<<"max_z = "<<max_z<<std::endl;
        std::cout<<"min_z = "<<min_z<<std::endl;
        std::cout<<"number of atoms in selection = "<<initPosSelection.atomCount()<<std::endl;

        // scale COM vector by total MASS:
//        centreOfMass[0] /= 1.0 * totalMass;
//        centreOfMass[1] /= 1.0 * totalMass;
//        centreOfMass[2] /= 1.0 * totalMass; 

        centreOfMass[0] = centreOfMass[0] / initPosSelection.atomCount();
        centreOfMass[1] = centreOfMass[1] / initPosSelection.atomCount();
        centreOfMass[2] = centreOfMass[2] / initPosSelection.atomCount(); 


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
 

        
   
    std::cout<<"vdwRadii.size() = "<<selVdwRadii.size()<<std::endl;
    std::cout<<"selection.atomCount() = "<<refSelection.atomCount()<<std::endl;


    pfm -> findPath();
 
    std::vector<gmx::RVec> path = pfm -> getPath();
    std::vector<real> radii = pfm -> getRadii();

    pfm.reset();
    std::cout<<"after pfm.reset()"<<std::endl;


   // const t_atoms *atoms;
   // snew(atoms);


    const char *outfile = "yo.pdb";
    const char *title = "title";
//    const t_atoms *atoms;
//    rvec tmp;
//    const rvec x[1] = {tmp};
//    const rvec *v;
//    int ePBC = 0;
//    const matrix box = NULL;

//    void write_sto_conf(outfile,
//                        title, 
//                        atoms, 
//                        x,
//                        v, 
//                        ePBC, 
//                        box);





    std::cout<<"path.size() = "<<path.size()<<std::endl; 
    // write path to DAT file:
    std::fstream datfile;
    std::string datfilename = "test.dat"; 
    datfile.open(datfilename.c_str(), std::fstream::out);
 
    for(unsigned int i=0; i<path.size(); i++)
    {
        datfile<<i<<"\t"                 // index
               <<path[i][0]<<"\t"        // x
               <<path[i][1]<<"\t"        // y
               <<path[i][2]<<"\t"        // z
               <<radii[i]                // radius
               <<std::endl;
    }


    datfile.close();


    // write path to PDB file:
    std::fstream pdbfile;
    std::string pdbfilename = "test.pdb"; 
    pdbfile.open(pdbfilename.c_str(), std::fstream::out);
    pdbfile.precision(3);

    pdbfile<<"HEADER"<<"   test"<<std::endl;
    pdbfile<<"TITLE"<<"   test"<<std::endl;
   

    for(unsigned int i=0; i<path.size(); i++)
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
            <<std::setw(8)<<path[i][0]*10.0        // x [Ang]
            <<std::setw(8)<<path[i][1]*10.0        // y [Ang]
            <<std::setw(8)<<path[i][2]*10.0        // z [Ang]
            <<std::setw(6)<<radii[i]*10.0          // occupancy [Ang]
            <<std::setw(6)<<radii[i]*10.0          // temperature factor [Ang]
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

    for(unsigned int i=0; i<path.size(); i++)
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


