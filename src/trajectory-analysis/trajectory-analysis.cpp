#include <algorithm>	// for std::max_element()
#include <cmath>		// for std::sqrt()
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <regex>

#include <gromacs/topology/atomprop.h>
#include <gromacs/random/threefry.h>
#include <gromacs/random/uniformrealdistribution.h>
#include <gromacs/fileio/confio.h>
#include <gromacs/utility/programcontext.h>

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include "trajectory-analysis/trajectory-analysis.hpp"

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_1D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"

#include "io/molecular_path_obj_exporter.hpp"
#include "io/json_doc_importer.hpp"
#include "io/analysis_data_json_exporter.hpp"

#include "trajectory-analysis/analysis_data_long_format_plot_module.hpp"
#include "trajectory-analysis/analysis_data_pdb_plot_module.hpp"

#include "path-finding/inplane_optimised_probe_path_finder.hpp"
#include "path-finding/optimised_direction_probe_path_finder.hpp"
#include "path-finding/naive_cylindrical_path_finder.hpp"
#include "path-finding/vdw_radius_provider.hpp"

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
    // HELP TEXT
    //-------------------------------------------------------------------------

	// set help text:
	static const char *const desc[] = {
		"This is a first prototype for the CHAP tool.",
		"There is NO HELP, you are on your own!"
	};
    settings -> setHelpText(desc);


    // SETTINGS
    //-------------------------------------------------------------------------

	// require the user to provide a topology file input:
    settings -> setFlag(TrajectoryAnalysisSettings::efRequireTop);

    // will not use periodic boundary conditions:
    settings -> setPBC(true);
    settings -> setFlag(TrajectoryAnalysisSettings::efNoUserPBC);

    // will make molecules whole:
    settings -> setRmPBC(false);
    settings -> setFlag(TrajectoryAnalysisSettings::efNoUserRmPBC);



    // OPTIONS
    //-------------------------------------------------------------------------

    // hardcoded defaults for multivalue options:
    std::vector<real> chanDirVec_ = {0.0, 0.0, 1.0};

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
    options -> addOption(RealOption("margin")
	                     .store(&poreMappingMargin_)
                         .defaultValue(1.5)
                         .description("Margin for residue mapping."));

    // file name for output json:
    options -> addOption(StringOption("json")
	                     .store(&jsonOutputFileName_)
                         .defaultValue("output.json")
                         .description("File name for JSON output."));


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
    options -> addOption(RealOption("pf-vdwr-fallback")
                         .store(&pfDefaultVdwRadius_)
                         .storeIsSet(&pfDefaultVdwRadiusIsSet_)
                         .defaultValue(-1.0)
                         .description("Fallback van-der-Waals radius for atoms that are not listed in van-der-Waals radius database"));
    const char * const allowedVdwRadiusDatabase[] = {"hole_amberuni",
                                                     "hole_bondi",
                                                     "hole_hardcore",
                                                     "hole_simple", 
                                                     "hole_xplor",
                                                     "user"};
    pfVdwRadiusDatabase_ = eVdwRadiusDatabaseHoleSimple;
    options -> addOption(EnumOption<eVdwRadiusDatabase>("pf-vdwr-database")
                         .enumValue(allowedVdwRadiusDatabase)
                         .store(&pfVdwRadiusDatabase_)
                         .description("Database of van-der-Waals radii to be used in pore finding"));
    options -> addOption(StringOption("pf-vdwr-json")
                         .store(&pfVdwRadiusJson_)
                         .storeIsSet(&pfVdwRadiusJsonIsSet_)
                         .description("User-defined set of van-der-Waals records in JSON format. Will be ignored unless -pf-vdwr-database is set to 'user'."));
    options -> addOption(StringOption("pf-method")
                         .store(&pfMethod_)
                         .defaultValue("inplane-optim")
                         .description("Path finding method. Only inplane-optim is implemented so far."));
    options -> addOption(RealOption("probe-step")
                         .store(&pfParams_["pfProbeStepLength"])
                         .defaultValue(0.025)
                         .description("Step length for probe movement. Defaults to 0.025 nm."));
    options -> addOption(RealOption("probe-radius")
                         .store(&pfParams_["pfProbeRadius"])
                         .defaultValue(0.0)
                         .description("Radius of probe. Defaults to 0.0, buggy for other values!"));
    options -> addOption(RealOption("max-free-dist")
                         .store(&pfParams_["pfProbeMaxRadius"])
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
    options -> addOption(Int64Option("sa-random-seed")
                         .store(&saRandomSeed_)
                         .storeIsSet(&saRandomSeedIsSet_)
                         .description("Seed for RNG used in simulated annealing. "));
    options -> addOption(IntegerOption("sa-max-cool")
                          .store(&saMaxCoolingIter_)
                          .defaultValue(1000)
                          .description("Maximum number of cooling iterations in one simulated annealing run. Defaults to 1000."));
    options -> addOption(IntegerOption("sa-cost-samples")
                         .store(&saNumCostSamples_)
                         .defaultValue(10)
                         .description("NOT IMPLEMENTED! Number of cost samples considered for convergence tolerance. Defaults to 10."));
    options -> addOption(RealOption("sa-conv-tol")
                         .store(&pfParams_["saConvTol"])
                         .defaultValue(1e-3)
                         .description("Relative tolerance for simulated annealing."));
    options -> addOption(RealOption("sa-init-temp")
                         .store(&pfParams_["saInitTemp"])
                         .defaultValue(0.1)
                         .description("Initital temperature for simulated annealing. Defaults to 0.1."));
    options -> addOption(RealOption("sa-cooling-fac")
                         .store(&pfParams_["saCoolingFactor"])
                         .defaultValue(0.98)
                         .description("Cooling factor using in simulated annealing. Defaults to 0.98."));
    options -> addOption(RealOption("sa-step")
                         .store(&pfParams_["saStepLengthFactor"])
                         .defaultValue(0.001)
                         .description("Step length factor used in candidate generation. Defaults to 0.001.")) ;
    options -> addOption(IntegerOption("nm-max-iter")
                         .store(&nmMaxIter_)
                         .defaultValue(100)
                         .description("Number of Nelder-Mead simplex iterations.")) ;
    options -> addOption(RealOption("nm-init-shift")
                         .store(&pfParams_["nmInitShift"])
                         .defaultValue(0.1)
                         .description("Distance of vertices in initial Nelder-Mead simplex.")) ;
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
    // PATH FINDING PARAMETERS
    //-------------------------------------------------------------------------

    // set inut-dependent defaults:
    if( !saRandomSeedIsSet_ )
    {
        saRandomSeed_ = gmx::makeRandomSeed();
    }


    std::cout<<"random seed = "<<saRandomSeed_<<std::endl;

    // set parameters in map:
    pfParams_["pfProbeMaxSteps"] = pfMaxProbeSteps_;

    pfParams_["pfCylRad"] = pfParams_["pfProbeMaxRadius"];
    pfParams_["pfCylNumSteps"] = pfParams_["pfProbeMaxSteps"];
    pfParams_["pfCylStepLength"] = pfParams_["pfProbeStepLength"];

    pfParams_["saMaxCoolingIter"] = saMaxCoolingIter_;
    pfParams_["saRandomSeed"] = saRandomSeed_;
    pfParams_["saNumCostSamples"] = saNumCostSamples_;

    pfParams_["nmMaxIter"] = nmMaxIter_;

	// set cutoff distance for grid search as specified in user input:
	nb_.setCutoff(cutoff_);
	std::cout<<"Setting cutoff to: "<<cutoff_<<std::endl;


    // PATH MAPPING PARAMETERS
    //-------------------------------------------------------------------------
    
    mappingParams_.nbhSearchCutoff_ = cutoff_ + poreMappingMargin_;
    mappingParams_.mapTol_ = 1e-7;
    mappingParams_.numPathSamples_ = 1000;
    mappingParams_.extrapDist_ = 100;


    // PREPARE DATSETS
    //-------------------------------------------------------------------------

    // multiple datasets:
    data_.setDataSetCount(4);
    DataSetNameList dataSetNames = {"path", "path.agg", "res.map", "solv.map"};
    ColumnHeaderList columnHeaders;

	// prepare data container for path data:
    data_.setColumnCount(0, 5);
    ColumnHeader columnHeaderPath = {"x", "y", "z", "s", "r"};
    columnHeaders.push_back(columnHeaderPath);

    // prepare data container for aggregate path data:
    data_.setColumnCount(1, 5);
    ColumnHeader columnHeaderAggregatePath = {"R.min", 
                                              "L", 
                                              "V", 
                                              "N",
                                              "N.sample"};
    columnHeaders.push_back(columnHeaderAggregatePath);

    // prepare data container for residue mapping:
    data_.setColumnCount(2, 9);
    ColumnHeader columnHeaderResMap = {"res.id", 
                                       "s", 
                                       "rho", 
                                       "phi", 
                                       "pl",
                                       "pf",
                                       "x",
                                       "y",
                                       "z"};
    columnHeaders.push_back(columnHeaderResMap);

    // prepare data container for solvent mapping:
    data_.setColumnCount(3, 9);
    ColumnHeader columnHeaderSolvMap = {"res.id", 
                                        "s", 
                                        "rho", 
                                        "phi", 
                                        "pore",
                                        "sample",
                                        "x",
                                        "y",
                                        "z"};
    columnHeaders.push_back(columnHeaderSolvMap);

    // prepare residue names:
    t_atoms allAtoms = top.topology() -> atoms;
    std::unordered_map<int, std::string> residueNames;
    std::cout<<"resinfo.nr = "<<allAtoms.nres<<std::endl;
    for(size_t i = 0; i < allAtoms.nres; i++)
    {
        residueNames[allAtoms.resinfo[i].nr] = std::string(*allAtoms.resinfo[i].name);        
    } 

    // add json exporter to data:
    AnalysisDataJsonExporterPointer jsonExporter(new AnalysisDataJsonExporter);
    jsonExporter -> setDataSetNames(dataSetNames);
    jsonExporter -> setColumnNames(columnHeaders);
    jsonExporter -> setResidueNames(residueNames);
    jsonExporter -> setFileName(jsonOutputFileName_);
    data_.addModule(jsonExporter);



    // RESIDUE MAPPING DATA
    //-------------------------------------------------------------------------


    // set dataset properties:
    dataResMapping_.setDataSetCount(1);
    dataResMapping_.setColumnCount(0, 6);   // refID s rho phi 
    dataResMapping_.setMultipoint(true);

    // add long format plot module:
    int j = 1;
    AnalysisDataLongFormatPlotModulePointer lfpltResMapping(new AnalysisDataLongFormatPlotModule(j));
    const char *fnResMapping = "res_mapping.dat";
    std::vector<char*> headerResMapping = {"t", 
                                           "mappedId", 
                                           "s", 
                                           "rho", 
                                           "phi", 
                                           "poreLining", 
                                           "poreFacing"};
    lfpltResMapping -> setFileName(fnResMapping);
    lfpltResMapping -> setHeader(headerResMapping);
    lfpltResMapping -> setPrecision(5);    // TODO: different treatment for integers?
    dataResMapping_.addModule(lfpltResMapping);


    // set pdb data set properties:
    dataResMappingPdb_.setDataSetCount(1);
    dataResMappingPdb_.setColumnCount(0, 7);

    // add PDB plot module:
    AnalysisDataPdbPlotModulePointer plotResMappingPdb(new AnalysisDataPdbPlotModule(3));
    plotResMappingPdb -> setFileName("res_mapping.pdb");
    dataResMappingPdb_.addModule(plotResMappingPdb);
    

    // PREPARE SELECTIONS FOR PORE PARTICLE MAPPING
    //-------------------------------------------------------------------------

    // prepare a centre of geometry selection collection:
    poreMappingSelCol_.setReferencePosType("res_cog");
    poreMappingSelCol_.setOutputPosType("res_cog");
  
    // selection of C-alpha atoms:
    // TODO: this will not work if only part of protein is specified as pore
    std::string refselSelText = refsel_.selectionText();
    std::string poreMappingSelCalString = "name CA";
    std::string poreMappingSelCogString = refselSelText;


    // create index groups from topology:
    // TODO: this will probably not work for custom index groups
    gmx_ana_indexgrps_t *poreIdxGroups;
    gmx_ana_indexgrps_init(&poreIdxGroups, 
                           top.topology(), 
                           NULL); 

    // create selections as defined above:
    poreMappingSelCal_ = poreMappingSelCol_.parseFromString(poreMappingSelCalString)[0];
    poreMappingSelCog_ = poreMappingSelCol_.parseFromString(poreMappingSelCogString)[0];
    poreMappingSelCol_.setTopology(top.topology(), 0);
    poreMappingSelCol_.setIndexGroups(poreIdxGroups);
    poreMappingSelCol_.compile();

    // free memory:
    gmx_ana_indexgrps_free(poreIdxGroups);

    // validate that there is a c-alpha for each residue:
    if( poreMappingSelCal_.posCount() != poreMappingSelCog_.posCount() )
    {
        std::cerr<<"ERROR: Could not find a C-alpha for each residue in pore forming group."
                 <<std::endl<<"Is your pore a protein?"<<std::endl;
        std::abort();
    }


    // PREPARE SELECTIONS FOR SOLVENT PARTICLE MAPPING
    //-------------------------------------------------------------------------

    // prepare centre of geometry seection collection:
    solvMappingSelCol_.setReferencePosType("res_cog");
    solvMappingSelCol_.setOutputPosType("res_cog");

    // create index groups from topology:
    // TODO: this will not work for custom index groups
    gmx_ana_indexgrps_t *solvIdxGroups;
    gmx_ana_indexgrps_init(&solvIdxGroups,
                           top.topology(),
                           NULL);

    // selection text:
    std::string solvMappingSelCogString = sel_[0].selectionText();

    // create selection as defined by user:
    solvMappingSelCog_ = solvMappingSelCol_.parseFromString(solvMappingSelCogString)[0];

    // compile the selections:
    solvMappingSelCol_.setTopology(top.topology(), 0);
    solvMappingSelCol_.setIndexGroups(solvIdxGroups);
    solvMappingSelCol_.compile();

    // free memory:
    gmx_ana_indexgrps_free(solvIdxGroups);

    
    // PREPARE TOPOLOGY QUERIES
    //-------------------------------------------------------------------------

	// load full topology:
	t_topology *topol = top.topology();	

	// access list of all atoms:
	t_atoms atoms = topol -> atoms;
    
	// create atomprop struct:
	gmx_atomprop_t aps = gmx_atomprop_init();

    
    // GET ATOM RADII FROM TOPOLOGY
    //-------------------------------------------------------------------------

    // get location of program binary from program context:
    const gmx::IProgramContext &programContext = gmx::getProgramContext();
    std::string radiusFilePath = programContext.fullBinaryPath();

    // obtain radius database location as relative path:
    auto lastSlash = radiusFilePath.find_last_of('/');
    radiusFilePath.replace(radiusFilePath.begin() + lastSlash + 1, 
                           radiusFilePath.end(), 
                           "data/vdwradii/");
        
    // select appropriate database file:
    if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseHoleAmberuni )
    {
        pfVdwRadiusJson_ = radiusFilePath + "hole_amberuni.json";
    }
    else if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseHoleBondi )
    {
        pfVdwRadiusJson_ = radiusFilePath + "hole_bondi.json";
    }
    else if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseHoleHardcore )
    {
        pfVdwRadiusJson_ = radiusFilePath + "hole_hardcore.json";
    }
    else if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseHoleSimple )
    {
        pfVdwRadiusJson_ = radiusFilePath + "hole_simple.json";
    }
    else if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseHoleXplor )
    {
        pfVdwRadiusJson_ = radiusFilePath + "hole_xplor.json";
    }
    else if( pfVdwRadiusDatabase_ == eVdwRadiusDatabaseUser )
    {
        // has user provided a file name?
        if( !pfVdwRadiusJsonIsSet_ )
        {
            std::cerr<<"ERROR: Option pfVdwRadiusDatabase set to 'user', but no custom van-der-Waals radii specified with pfVdwRadiusJson."<<std::endl;
            std::abort();
        }
    }

    // import vdW radii JSON: 
    JsonDocImporter jdi;
    rapidjson::Document radiiDoc = jdi(pfVdwRadiusJson_.c_str());
   
    // create radius provider and build lookup table:
    VdwRadiusProvider vrp;
    try
    {
        vrp.lookupTableFromJson(radiiDoc);
    }
    catch( std::exception& e )
    {
        std::cerr<<"ERROR while creating van der Waals radius lookup table:"<<std::endl;
        std::cerr<<e.what()<<std::endl; 
        std::abort();
    }

	// find largest van der Waals radius in system:
//	maxVdwRadius_ = *std::max_element(vdwRadii_.begin(), vdwRadii_.end());


    // TRACK C-ALPHAS and RESIDUE INDICES
    //-------------------------------------------------------------------------
   
    // loop through all atoms, get index lists for c-alphas and residues:
    for(int i = 0; i < atoms.nr; i++)
    {
        // check for calpha:
        if( std::strcmp(*atoms.atomname[i], "CA") == 0 )
        {
           poreCAlphaIndices_.push_back(i); 
        }

        // track residue ID of atoms: 
        residueIndices_.push_back(atoms.atom[i].resind);
        atomResidueMapping_[i] = atoms.atom[i].resind;
        residueAtomMapping_[atoms.atom[i].resind].push_back(i);
    }

    // remove duplicate residue indices:
    std::vector<int>::iterator it;
    it = std::unique(residueIndices_.begin(), residueIndices_.end());
    residueIndices_.resize(std::distance(residueIndices_.begin(), it));

    //
    ConstArrayRef<int> refselAtomIdx = refsel_.atomIndices();
    for(it = residueIndices_.begin(); it != residueIndices_.end(); it++)
    {
        // current residue id:
        int resId = *it;
    
        // get vector of all atom indices in this residue:
        std::vector<int> atomIdx = residueAtomMapping_[resId];

        // for each atom in residue, check if it belongs to pore selection:
        bool addResidue = false;
        std::vector<int>::iterator jt;
        for(jt = atomIdx.begin(); jt != atomIdx.end(); jt++)
        {
            // check if atom belongs to pore selection:
            if( std::find(refselAtomIdx.begin(), refselAtomIdx.end(), *jt) != refselAtomIdx.end() )
            {
                // add atom to list of pore atoms:
                poreAtomIndices_.push_back(*jt);

                // if at least one atom belongs to pore group, the whole residue will be considered:
                addResidue = true;
            }
        }

        // add residue, if at least one atom is part of pore:
        if( addResidue == true )
        {
            poreResidueIndices_.push_back(resId);
        }
    }


    // FINALISE ATOMPROP QUERIES
    //-------------------------------------------------------------------------
    
	// delete atomprop struct:
	gmx_atomprop_destroy(aps);

    // set user-defined default radius?
    if( pfDefaultVdwRadiusIsSet_ )
    {
        vrp.setDefaultVdwRadius(pfDefaultVdwRadius_);
    }

    // build vdw radius lookup map:
    try
    {
        vdwRadii_ = vrp.vdwRadiiForTopology(top, refsel_.mappedIds());
    }
    catch( std::exception& e )
    {
        std::cerr<<"ERROR in van der Waals radius lookup:"<<std::endl;
        std::cerr<<e.what()<<std::endl;
        std::abort();
    } 

    // find maximum van der Waals radius:
    maxVdwRadius_ = std::max_element(vdwRadii_.begin(), vdwRadii_.end()) -> second;
}




/*
 *
 */
void
trajectoryAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                 TrajectoryAnalysisModuleData *pdata)
{
	// get thread-local selections:
	const Selection &refSelection = pdata -> parallelSelection(refsel_);
    const Selection &initProbePosSelection = pdata -> parallelSelection(initProbePosSelection_);

    // get data handles for this frame:
	AnalysisDataHandle dh = pdata -> dataHandle(data_);
    AnalysisDataHandle dhResMapping = pdata -> dataHandle(dataResMapping_);

	// get data for frame number frnr into data handle:
    dh.startFrame(frnr, fr.time);
    dhResMapping.startFrame(frnr, fr.time);


    // UPDATE INITIAL PROBE POSITION FOR THIS FRAME
    //-------------------------------------------------------------------------

    // recalculate initial probe position based on reference group COG:
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


    // GET VDW RADII FOR SELECTION
    //-------------------------------------------------------------------------
    // TODO: Move this to separate class and test!
    // TODO: Should then also work for coarse-grained situations!

	// create vector of van der Waals radii and allocate memory:
    std::vector<real> selVdwRadii;
	selVdwRadii.reserve(refSelection.atomCount());

    // loop over all atoms in system and get vdW-radii:
	for(int i=0; i<refSelection.atomCount(); i++)
    {
        // get global index of i-th atom in selection:
        gmx::SelectionPosition atom = refSelection.position(i);
        int idx = atom.mappedId();

		// add radius to vector of radii:
		selVdwRadii.push_back(vdwRadii_.at(idx));
	}


	// PORE FINDING AND RADIUS CALCULATION
	// ------------------------------------------------------------------------

    // vectors as RVec:
    RVec initProbePos(pfInitProbePos_[0], pfInitProbePos_[1], pfInitProbePos_[2]);
    RVec chanDirVec(pfChanDirVec_[0], pfChanDirVec_[1], pfChanDirVec_[2]); 

    // create path finding module:
    std::unique_ptr<AbstractPathFinder> pfm;
    if( pfMethod_ == "inplane-optim" )
    {
        // create inplane-optimised path finder:
        pfm.reset(new InplaneOptimisedProbePathFinder(pfParams_,
                                                      initProbePos,
                                                      chanDirVec,
                                                      *pbc,
                                                      refSelection,
                                                      selVdwRadii));
    }
    else if( pfMethod_ == "optim-direction" )
    {
        std::cerr<<"ERROR: Optimised direction path finding is not implemented!"<<std::endl;
        std::abort();
    }   
    else if( pfMethod_ == "naive-cylindrical" )
    {        
        // create the naive cylindrical path finder:
        pfm.reset(new NaiveCylindricalPathFinder(pfParams_,
                                                 initProbePos,
                                                 chanDirVec));
    }











    std::cout<<std::endl;
    std::cout<<"initProbePos ="<<" "
             <<pfInitProbePos_[0]<<" "
             <<pfInitProbePos_[1]<<" "
             <<pfInitProbePos_[2]<<" "
             <<std::endl;







    // PATH FINDING
    //-------------------------------------------------------------------------


    // run path finding algorithm on current frame:
    std::cout<<"finding permeation pathway ... ";
    std::cout.flush();
    clock_t tPathFinding = std::clock();
    pfm -> findPath();
    tPathFinding = (std::clock() - tPathFinding)/CLOCKS_PER_SEC;
    std::cout<<"done in  "<<tPathFinding<<" sec"<<std::endl;


    // retrieve molecular path object:
    std::cout<<"preparing pathway object ... ";
    clock_t tMolPath = std::clock();
    MolecularPath molPath = pfm -> getMolecularPath();
    tMolPath = (std::clock() - tMolPath)/CLOCKS_PER_SEC;
    std::cout<<"done in  "<<tMolPath<<" sec"<<std::endl;
    
    std::vector<gmx::RVec> pathPoints = molPath.pathPoints();
    std::vector<real> pathRadii = molPath.pathRadii();

    // add path data to data handle:
    dh.selectDataSet(0);

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



    // MAP PORE PARTICLES ONTO PATHWAY
    //-------------------------------------------------------------------------

    std::cout<<std::endl;

 
    // evaluate pore mapping selection for this frame:
    t_trxframe frame = fr;
    poreMappingSelCol_.evaluate(&frame, pbc);
    const gmx::Selection poreMappingSelCal = pdata -> parallelSelection(poreMappingSelCal_);    
    const gmx::Selection poreMappingSelCog = pdata -> parallelSelection(poreMappingSelCog_);    


    // map pore residue COG onto pathway:
    std::cout<<"mapping pore residue COG onto pathway ... ";
    clock_t tMapResCog = std::clock();
    std::map<int, gmx::RVec> poreCogMappedCoords = molPath.mapSelection(
            poreMappingSelCog, 
            mappingParams_,
            pbc);
    tMapResCog = (std::clock() - tMapResCog)/CLOCKS_PER_SEC;
    std::cout<<"mapped "<<poreCogMappedCoords.size()
             <<" particles in "<<1000*tMapResCog<<" ms"<<std::endl;

    // map pore residue C-alpha onto pathway:
    std::cout<<"mapping pore residue C-alpha onto pathway ... ";
    clock_t tMapResCal = std::clock();
    std::map<int, gmx::RVec> poreCalMappedCoords = molPath.mapSelection(
            poreMappingSelCal, 
            mappingParams_,
            pbc);
    tMapResCal = (std::clock() - tMapResCal)/CLOCKS_PER_SEC;
    std::cout<<"mapped "<<poreCalMappedCoords.size()
             <<" particles in "<<1000*tMapResCal<<" ms"<<std::endl;

    
    // check if particles are pore-lining:
    std::cout<<"checking which residues are pore-lining ... ";
    clock_t tResPoreLining = std::clock();
    std::map<int, bool> poreLining = molPath.checkIfInside(poreCogMappedCoords, poreMappingMargin_);
    int nPoreLining = 0;
    for(auto jt = poreLining.begin(); jt != poreLining.end(); jt++)
    {
        if( jt -> second == true )
        {
            nPoreLining++;
        }
    }
    tResPoreLining = (std::clock() - tResPoreLining)/CLOCKS_PER_SEC;
    std::cout<<"found "<<nPoreLining<<" pore lining residues in "
             <<1000*tResPoreLining<<" ms"<<std::endl;

    // check if residues are pore-facing:
    // TODO: make this conditional on whether C-alphas are available
    
    std::cout<<"checking which residues are pore-facing ... ";
    clock_t tResPoreFacing = std::clock();
    std::map<int, bool> poreFacing;
    int nPoreFacing = 0;
    for(auto it = poreCogMappedCoords.begin(); it != poreCogMappedCoords.end(); it++)
    {
        // is residue pore lining and has COG closer to centreline than CA?
        if( it -> second[1] > poreCalMappedCoords[it->first][1] &&
            poreLining[it -> first] == true )
        {
            poreFacing[it->first] = true;
            nPoreFacing++;
        }
        else
        {
            poreFacing[it->first] = false;            
        }
    }
    tResPoreFacing = (std::clock() - tResPoreFacing)/CLOCKS_PER_SEC;
    std::cout<<"found "<<nPoreFacing<<" pore facing residues in "
             <<1000*tResPoreFacing<<" ms"<<std::endl;
    

    // add points inside to data frame:
    for(auto it = poreCogMappedCoords.begin(); it != poreCogMappedCoords.end(); it++)
    {
        SelectionPosition pos = poreMappingSelCog.position(it->first);
        
        // add points to dataset:
        dhResMapping.setPoint(0, pos.mappedId());                    // refId
        dhResMapping.setPoint(1, it -> second[0]); // s
        dhResMapping.setPoint(2, it -> second[1]); // rho
        dhResMapping.setPoint(3, it -> second[2]);// phi
        dhResMapping.setPoint(4, poreLining[it -> first]);             // poreLining
        dhResMapping.setPoint(5, poreFacing[it -> first]);             // poreFacing TODO
        dhResMapping.finishPointSet();
    }
    
    // now add mapped residue coordinates to data handle:
    dh.selectDataSet(2);
    
    // add mapped residues to data container:
    for(auto it = poreCogMappedCoords.begin(); it != poreCogMappedCoords.end(); it++)
    {
         dh.setPoint(0, poreMappingSelCog.position(it -> first).mappedId()); // res.id
         dh.setPoint(1, it -> second[0]);     // s
         dh.setPoint(2, it -> second[1]);     // rho
         dh.setPoint(3, it -> second[3]);     // phi
         dh.setPoint(4, poreLining[it -> first]);     // pore lining?
         dh.setPoint(5, poreFacing[it -> first]);     // pore facing?
         dh.setPoint(6, poreMappingSelCog.position(it -> first).x()[0]);  // x
         dh.setPoint(7, poreMappingSelCog.position(it -> first).x()[1]);  // y
         dh.setPoint(8, poreMappingSelCog.position(it -> first).x()[2]);  // z
         dh.finishPointSet();
    }


    // MAP SOLVENT PARTICLES ONTO PATHWAY
    //-------------------------------------------------------------------------

    // evaluate solevnt mapping selections for this frame:
    t_trxframe tmpFrame = fr;
    solvMappingSelCol_.evaluate(&tmpFrame, pbc);

    // TODO: make this a parameter:
    real solvMappingMargin_ = 0.0;
        
    // get thread-local selection data:
    const Selection solvMapSel = pdata -> parallelSelection(solvMappingSelCog_);

    // map particles onto pathway:
    std::cout<<"mapping solvent particles onto pathway ... ";
    clock_t tMapSol = std::clock();
    std::map<int, gmx::RVec> solventMappedCoords = molPath.mapSelection(
            solvMapSel, 
            mappingParams_,
            pbc);
    tMapSol = (std::clock() - tMapSol)/CLOCKS_PER_SEC;
    std::cout<<"mapped "<<solventMappedCoords.size()
             <<" particles in "<<1000*tMapSol<<" ms"<<std::endl;


    std::cout<<"solvMapSel.posCount = "<<solvMapSel.posCount()<<std::endl;
    std::cout<<"solventMappedCoords = "<<solventMappedCoords.size()<<std::endl;


    // find particles inside path (i.e. pore plus bulk sampling regime):
    std::cout<<"finding solvent particles inside pore ... ";
    clock_t tSolInsideSample = std::clock();
    std::map<int, bool> solvInsideSample = molPath.checkIfInside(
            solventMappedCoords, solvMappingMargin_);
    int numSolvInsideSample = 0;
    for(auto jt = solvInsideSample.begin(); jt != solvInsideSample.end(); jt++)
    {            
        if( jt -> second == true )
        {
            numSolvInsideSample++;
        }
    }
    tSolInsideSample = (std::clock() - tSolInsideSample)/CLOCKS_PER_SEC;
    std::cout<<"found "<<numSolvInsideSample<<" solvent particles inside pore in "
             <<1000*tSolInsideSample<<" ms"<<std::endl;


    // find particles inside pore:
    std::cout<<"finding solvent particles inside pore ... ";
    clock_t tSolInsidePore = std::clock();
    std::map<int, bool> solvInsidePore = molPath.checkIfInside(
            solventMappedCoords, 
            solvMappingMargin_,
            molPath.sLo(),
            molPath.sHi());
    int numSolvInsidePore = 0;
    for(auto jt = solvInsidePore.begin(); jt != solvInsidePore.end(); jt++)
    {            
        if( jt -> second == true )
        {
            numSolvInsidePore++;
        }
    }
    tSolInsidePore = (std::clock() - tSolInsidePore)/CLOCKS_PER_SEC;
    std::cout<<"found "<<numSolvInsidePore<<" solvent particles inside pore in "
             <<1000*tSolInsidePore<<" ms"<<std::endl;


    // now add mapped residue coordinates to data handle:
    dh.selectDataSet(3);
    
    // add mapped residues to data container:
    for(auto it = solventMappedCoords.begin(); 
        it != solventMappedCoords.end(); 
        it++)
    {
         dh.setPoint(0, solvMapSel.position(it -> first).mappedId()); // res.id
         dh.setPoint(1, it -> second[0]);     // s
         dh.setPoint(2, it -> second[1]);     // rho
         dh.setPoint(3, it -> second[3]);     // phi
         dh.setPoint(4, solvInsidePore[it -> first]);     // inside pore
         dh.setPoint(5, solvInsideSample[it -> first]);     // inside sample
         dh.setPoint(6, solvMapSel.position(it -> first).x()[0]);  // x
         dh.setPoint(7, solvMapSel.position(it -> first).x()[1]);  // y
         dh.setPoint(8, solvMapSel.position(it -> first).x()[2]);  // z
         dh.finishPointSet();
    }


    // ADD AGGREGATE DATA TO PARALLELISABLE CONTAINER
    //-------------------------------------------------------------------------


    // add aggegate path data:
    dh.selectDataSet(1);

    // only one point per frame:
    dh.setPoint(0, molPath.minRadius().second);
    dh.setPoint(1, molPath.length());
    dh.setPoint(2, molPath.volume());
    dh.setPoint(3, numSolvInsidePore); 
    dh.setPoint(4, numSolvInsideSample); 
    dh.finishPointSet();
    


    // WRITE PORE TO OBJ FILE
    //-------------------------------------------------------------------------

    MolecularPathObjExporter molPathExp;
    molPathExp("pore.obj",
               molPath);


    // FINISH FRAME
    //-------------------------------------------------------------------------

    std::cout<<std::endl;

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
    std::cout<<"finished analysis"<<std::endl;
}




void
trajectoryAnalysis::writeOutput()
{
    std::cout<<"datSetCount = "<<data_.dataSetCount()<<std::endl
             <<"columnCount = "<<data_.columnCount()<<std::endl
             <<"frameCount = "<<data_.frameCount()<<std::endl
             <<std::endl;
}


