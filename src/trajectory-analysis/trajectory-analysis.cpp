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
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"

#include "trajectory-analysis/trajectory-analysis.hpp"

#include "aggregation/number_density_calculator.hpp"
#include "aggregation/boltzmann_energy_calculator.hpp"

#include "config/version.hpp"

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_1D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"

#include "io/molecular_path_obj_exporter.hpp"
#include "io/json_doc_importer.hpp"
#include "io/analysis_data_json_frame_exporter.hpp"
#include "io/summary_statistics_json_converter.hpp"

#include "statistics/summary_statistics.hpp"
#include "statistics/histogram_density_estimator.hpp"
#include "statistics/kernel_density_estimator.hpp"

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
    , pfProbeRadius_(0.0)
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
    // register dataset:
    // TODO: this data set should be made obsolete
    registerAnalysisDataset(&dataResMapping_, "resMapping");


    registerAnalysisDataset(&frameStreamData_, "frameStreamData");
    frameStreamData_.setMultipoint(true); 



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



    // GENERAL OPTIONS
    //-------------------------------------------------------------------------

	options -> addOption(SelectionOption("sel-pathway")
	                     .store(&refsel_).required()
		                 .description("Reference group that defines the "
                                      "permeation pathway (usually "
                                      "'Protein') "));

	options -> addOption(SelectionOption("sel-solvent")
                         .storeVector(&sel_).required()
	                     .description("Group of small particles to calculate "
                                      "density of (usually 'Water')"));


    // OUTPUT OPTIONS
    // ------------------------------------------------------------------------

    // TODO: remove
    options -> addOption(StringOption("ppfn")
                         .store(&poreParticleFileName_)
                         .defaultValue("pore_particles.dat")
                         .description("Name of file containing pore particle "
                                      "positions over time."));

    // TODO remove
    options -> addOption(StringOption("spfn")
                         .store(&smallParticleFileName_)
                         .defaultValue("small_particles.dat")
                         .description("Name of file containing small particle "
                                      "positions (i.e. water particle "
                                      "positions) over time."));

    // TODO: find better solution for this
    options -> addOption(StringOption("o")
                         .store(&poreProfileFileName_)
                         .defaultValue("pore_profile.dat")
                         .description("Name of file containing pore radius, "
                                      "small particle density, and small "
                                      "particle energy as function of the "
                                      "permeation coordinate."));

    // TODO: is this used?
    options -> addOption(IntegerOption("num-out-pts")
                         .store(&nOutPoints_)
                         .defaultValue(1000)
                         .description("Number of sample points of pore centre "
                                      "line that are written to output."));

    // TODO: more exressive or shorter command line flag?
    options -> addOption(StringOption("json")
	                     .store(&jsonOutputFileName_)
                         .defaultValue("output.json")
                         .description("File name for JSON output."));

    // TODO find better solution for this
    options -> addOption(StringOption("obj")
	                     .store(&objOutputFileName_)
                         .defaultValue("output.obj")
                         .description("File name for OBJ output (testing)."));

    // TODO find better solution for this.
    options -> addOption(StringOption("resm")
	                     .store(&resMappingOutFileName_)
                         .defaultValue("res_mapping.dat")
                         .description("Residue mapping data (testing)."));


    // PATH FINDING PARAMETERS
    //-------------------------------------------------------------------------

    // TODO: make this enum
    options -> addOption(StringOption("pf-method")
                         .store(&pfMethod_)
                         .defaultValue("inplane-optim")
                         .description("Path finding method. Only "
                                      "inplane-optim is implemented so far."));

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
                         .description("Database of van-der-Waals radii to be "
                                      "used in pore finding"));

    options -> addOption(RealOption("pf-vdwr-fallback")
                         .store(&pfDefaultVdwRadius_)
                         .storeIsSet(&pfDefaultVdwRadiusIsSet_)
                         .defaultValue(-1.0)
                         .description("Fallback van-der-Waals radius for "
                                      "atoms that are not listed in "
                                      "van-der-Waals radius database. If "
                                      "negative, an error will be thrown if "
                                      "the database does not contain a "
                                      "van-der-Waals radii for all particles "
                                      "in the pathway defining group."));

    options -> addOption(StringOption("pf-vdwr-json")
                         .store(&pfVdwRadiusJson_)
                         .storeIsSet(&pfVdwRadiusJsonIsSet_)
                         .description("JSON file with user defined "
                                      "van-der-Waals radii. Will be "
                                      "ignored unless -pf-vdwr-database is "
                                      "set to 'user'."));

    const char * const allowedPathAlignmentMethod[] = {"none",
                                                       "ipp"};
    pfPathAlignmentMethod_ = ePathAlignmentMethodIpp;
    options -> addOption(EnumOption<ePathAlignmentMethod>("pf-align-method")
                         .enumValue(allowedPathAlignmentMethod)
                         .store(&pfPathAlignmentMethod_)
                         .description("Method for aligning pathway "
                                      "coordinates across time steps"));

    options -> addOption(RealOption("pf-probe-step")
                         .store(&pfProbeStepLength_)
                         .defaultValue(0.025)
                         .description("Step length for probe movement."));

    // TODO: remove this
    options -> addOption(RealOption("pf-probe-radius")
                         .store(&pfPar_["pfProbeRadius"])
                         .defaultValue(0.0)
                         .description("BUGGY! Radius of probe."));

    options -> addOption(RealOption("pf-max-free-dist")
                         .store(&pfMaxProbeRadius_)
                         .defaultValue(1.0)
                         .description("Maximum radius of pore."));

    options -> addOption(IntegerOption("pf-max-probe-steps")
                         .store(&pfMaxProbeSteps_)
                         .defaultValue(10000)
                         .description("Maximum number of steps the probe is "
                                      "moved in either direction."));

    options -> addOption(SelectionOption("pf-sel-ipp")
                         .store(&ippsel_)
                         .storeIsSet(&ippselIsSet_)
	                     .description("Reference group from which to "
                                      "determine the initial probe position "
                                      "for the path finding algorithm. If "
                                      "unspecified, this defaults to the "
                                      "overall path defining group. Will be "
                                      "overridden if init-probe-pos is set "
                                      "explicitly."));

    options -> addOption(RealOption("pf-init-probe-pos")
                         .storeVector(&pfInitProbePos_)
                         .storeIsSet(&pfInitProbePosIsSet_)
                         .valueCount(3)
                         .description("Initial position of probe in "
                                      "probe-based pore finding algorithms. "
                                      "If set explicitly, it will overwrite "
                                      "the COM-based initial position set "
                                      "with the ippselflag."));

    std::vector<real> chanDirVec_ = {0.0, 0.0, 1.0};
    options -> addOption(RealOption("pf-chan-dir-vec")
                         .storeVector(&pfChanDirVec_)
                         .storeIsSet(&pfChanDirVecIsSet_)
                         .valueCount(3)
                         .description("Channel direction vector. Will be "
                                      "normalised to unit vector internally. "
                                      "If unset pore is assumed to be "
                                      "oriented in z-direction."));
   
    // TODO: should be possible to determine this automatically from
    // max-free-dist and largest vdW radius
    options -> addOption(DoubleOption("pf-cutoff")
	                     .store(&cutoff_)
                         .storeIsSet(&cutoffIsSet_)
                         .defaultValue(0.0)
                         .description("Cutoff for distance searches in path "
                                      "finding algorithm. A value of zero "
                                      "or less means no cutoff is applied."));
 


    // OPTIMISATION PARAMETERS
    //-------------------------------------------------------------------------

    options -> addOption(Int64Option("sa-seed")
                         .store(&saRandomSeed_)
                         .storeIsSet(&saRandomSeedIsSet_)
                         .description("Seed used in pseudo random number "
                                      "generation for simulated annealing. "
                                      "If not set explicitly, a random seed "
                                      "is used."));

    options -> addOption(IntegerOption("sa-max-iter")
                          .store(&saMaxCoolingIter_)
                          .defaultValue(0)
                          .description("Maximum number of cooling iterations "
                                       "in one simulated annealing run."));
                          
    // TODO remove this
    options -> addOption(IntegerOption("sa-cost-samples")
                         .store(&saNumCostSamples_)
                         .defaultValue(10)
                         .description("NOT IMPLEMENTED! Number of cost samples"
                                      " considered for convergence tolerance."));

    // TODO remove this
    options -> addOption(RealOption("sa-conv-tol")
                         .store(&pfPar_["saConvTol"])
                         .defaultValue(1e-3)
                         .description("Simulated annealing relative tolerance."));

    options -> addOption(RealOption("sa-init-temp")
                         .store(&pfPar_["saInitTemp"])
                         .defaultValue(0.1)
                         .description("Simulated annealing initial "
                                      "temperature."));

    options -> addOption(RealOption("sa-cooling-fac")
                         .store(&pfPar_["saCoolingFactor"])
                         .defaultValue(0.98)
                         .description("Simulated annealing cooling factor."));

    options -> addOption(RealOption("sa-step")
                         .store(&pfPar_["saStepLengthFactor"])
                         .defaultValue(0.001)
                         .description("Step length factor used in candidate "
                                      "generation. Defaults to 0.001."));

    options -> addOption(IntegerOption("nm-max-iter")
                         .store(&nmMaxIter_)
                         .defaultValue(100)
                         .description("Number of Nelder-Mead simplex "
                                      "iterations in path finding algorithm."));

    options -> addOption(RealOption("nm-init-shift")
                         .store(&pfPar_["nmInitShift"])
                         .defaultValue(0.1)
                         .description("Distance of vertices in initial "
                                      "Nelder-Mead simplex."));


    // PATH MAPPING PARAMETERS
    //-------------------------------------------------------------------------

    options -> addOption(RealOption("pm-tol")
                         .store(&mappingParams_.mapTol_)
                         .defaultValue(1e-7)
                         .description("Tolerance threshold for mapping "
                                      "particles onto molecular pathway."));

    options -> addOption(RealOption("pm-extrap-dist")
                         .store(&mappingParams_.extrapDist_)
                         .defaultValue(10)
                         .description("Extrapolation distance for sampling "
                                      "path points outside the pore when "
                                      "mapping particles onto molecular "
                                      "pathway. Should be large enough to "
                                      "sample a sufficient portion of the "
                                      "bulk regime to reliably estimate bulk "
                                      "solvent density."));

    options -> addOption(RealOption("pm-sample-step")
                         .store(&mappingParams_.sampleStep_)
                         .defaultValue(0.01)
                         .description("Arc length distance of path samples "
                                      "when mapping particles onto molecular "
                                      "pathway."));

    options -> addOption(RealOption("pm-pl-margin")
	                     .store(&poreMappingMargin_)
                         .defaultValue(0.9)
                         .description("Margin for determining pathway lining "
                                      "residues. A residue is considered to "
                                      "be pathway lining if it is no further "
                                      "than the local path radius plus this "
                                      "margin from the pathway's centre "
                                      "line."));


    // DENSITY ESTIMATION PARAMETERS
    //-------------------------------------------------------------------------

    const char * const allowedDensityEstimationMethod[] = {"histogram",
                                                           "kernel"};
    deMethod_ = eDensityEstimatorKernel;
    options -> addOption(EnumOption<eDensityEstimator>("de-method")
                         .enumValue(allowedDensityEstimationMethod)
                         .store(&deMethod_)
                         .description("Method used for estimating the "
                                      "probability density of the solvent "
                                      "particles along the permeation "
                                      "pathway"));
    
    options -> addOption(RealOption("de-res")
                         .store(&deResolution_)
                         .defaultValue(0.01)
                         .description("Spatial resolution of the density "
                                      "estimator. In case of a historgam, "
                                      "this is the bin widtj, in case of a "
                                      "kernel density estimator, this is the "
                                      "spacing of the evaluation points."));

    // TODO add functionality to determine this automatically
    options -> addOption(RealOption("de-bandwidth")
                         .store(&deBandWidth_)
                         .defaultValue(0.1)
                         .description("Bandwidth for the kernel density "
                                      "estimator. Ignored for other "
                                      "methods."));

    options -> addOption(RealOption("de-eval-cutoff")
                         .store(&deEvalRangeCutoff_)
                         .defaultValue(5)
                         .description("Evaluation range cutoff for kernel "
                                      "density estimator in multiples of "
                                      "bandwidth. Ignored for other methods. "
                                      "Ensures that the density falls off "
                                      "smoothly to zero outside the data "
                                      "range."));


    // MISC PARAMETERS
    //-------------------------------------------------------------------------

    // TODO: remove this
    options -> addOption(BooleanOption("debug-output")
                         .store(&debug_output_)
                         .description("When this flag is set, the program "
                                      "will write additional information."));
}




/*
 * 
 */
void
trajectoryAnalysis::initAnalysis(const TrajectoryAnalysisSettings& /*settings*/,
                                 const TopologyInformation &top)
{
    // PATH FINDING PARAMETERS
    //-------------------------------------------------------------------------

    // set inut-dependent defaults:
    if( !saRandomSeedIsSet_ )
    {
        saRandomSeed_ = gmx::makeRandomSeed();
    }

    // set parameters in map:
    pfPar_["pfProbeMaxSteps"] = pfMaxProbeSteps_;

    pfPar_["pfCylRad"] = pfMaxProbeRadius_;
    pfPar_["pfCylNumSteps"] = pfPar_["pfProbeMaxSteps"];
    pfPar_["pfCylStepLength"] = pfProbeStepLength_;

    pfPar_["saMaxCoolingIter"] = saMaxCoolingIter_;
    pfPar_["saRandomSeed"] = saRandomSeed_;
    pfPar_["saNumCostSamples"] = saNumCostSamples_;

    pfPar_["nmMaxIter"] = nmMaxIter_;


    // 
    pfParams_.setProbeStepLength(pfProbeStepLength_);
    pfParams_.setMaxProbeRadius(pfMaxProbeRadius_);
    pfParams_.setMaxProbeSteps(pfMaxProbeSteps_);
    
    if( cutoffIsSet_ )
    {
        pfParams_.setNbhCutoff(cutoff_);
    }

    // PATH MAPPING PARAMETERS
    //-------------------------------------------------------------------------

    // sanity checks and automatic defaults:
    if( mappingParams_.mapTol_ <= 0.0 )
    {
        throw(std::runtime_error("Mapping tolerance parameter pm-tol must be positive."));
    }

    if( mappingParams_.extrapDist_ <= 0 )
    {
        throw(std::runtime_error("Extrapolation distance set with pm-extrap-dist may not be negative."));
    }

    if( mappingParams_.sampleStep_ <= 0 )
    {
        throw(std::runtime_error("Sampling step set with pm-sample-step must be positive."));
    }


    // DENSITY ESTIMATION PARAMETERS
    //-------------------------------------------------------------------------

    // which estimator will be used?
    if( deMethod_ == eDensityEstimatorHistogram )
    {
        deParams_.setBinWidth(deResolution_);
    }
    else if( deMethod_ == eDensityEstimatorKernel )
    {
        deParams_.setKernelFunction(eKernelFunctionGaussian);
        deParams_.setBandWidth(deBandWidth_);
        deParams_.setEvalRangeCutoff(deEvalRangeCutoff_);
        deParams_.setMaxEvalPointDist(deResolution_);
    }


    // PREPARE DATSETS
    //-------------------------------------------------------------------------

    // prepare per frame data stream:
    frameStreamData_.setDataSetCount(7);
    std::vector<std::string> frameStreamDataSetNames = {
            "pathSummary",
            "molPathOrigPoints",
            "molPathRadiusSpline",
            "molPathCentreLineSpline",
            "residuePositions",
            "solventPositions",
            "solventDensitySpline"};
    std::vector<std::vector<std::string>> frameStreamColumnNames;


    // prepare container for aggregated data:
    frameStreamData_.setColumnCount(0, 7);
    frameStreamColumnNames.push_back({"minRadius",
                                      "length",
                                      "volume",
                                      "numPath",
                                      "numSample",
                                      "solventRangeLo",
                                      "solventRangeHi"});

    // prepare container for original path points:
    frameStreamData_.setColumnCount(1, 4);
    frameStreamColumnNames.push_back({"x", 
                                      "y",
                                      "z",
                                      "r"});

    // prepare container for path radius:
    frameStreamData_.setColumnCount(2, 2);
    frameStreamColumnNames.push_back({"knots", 
                                      "ctrl"});

    // prepare container for pathway spline:
    frameStreamData_.setColumnCount(3, 4);
    frameStreamColumnNames.push_back({"knots", 
                                      "ctrlX",
                                      "ctrlY",
                                      "ctrlZ"});

    // prepare container for residue mapping results:
    frameStreamData_.setColumnCount(4, 9);
    frameStreamColumnNames.push_back({"resId",
                                      "s",
                                      "rho",
                                      "phi",
                                      "poreLining",
                                      "poreFacing",
                                      "x",
                                      "y",
                                      "z"});

    // prepare container for solvent mapping:
    frameStreamData_.setColumnCount(5, 9);
    frameStreamColumnNames.push_back({"resId", 
                                      "s",
                                      "rho",
                                      "phi",
                                      "inPore",
                                      "inSample",
                                      "x",
                                      "y",
                                      "z"});

    // prepare container for solvent density:
    frameStreamData_.setColumnCount(6, 2);
    frameStreamColumnNames.push_back({"knots", 
                                      "ctrl"});

    // add JSON exporter to frame stream data:
    AnalysisDataJsonFrameExporterPointer jsonFrameExporter(new AnalysisDataJsonFrameExporter);
    jsonFrameExporter -> setDataSetNames(frameStreamDataSetNames);
    jsonFrameExporter -> setColumnNames(frameStreamColumnNames);
    std::string frameStreamFileName = std::string("stream_") + jsonOutputFileName_;
    jsonFrameExporter -> setFileName(frameStreamFileName);
    frameStreamData_.addModule(jsonFrameExporter);


    // RESIDUE MAPPING DATA
    //-------------------------------------------------------------------------


    // set dataset properties:
    dataResMapping_.setDataSetCount(1);
    dataResMapping_.setColumnCount(0, 6);   // refID s rho phi 
    dataResMapping_.setMultipoint(true);

    // add long format plot module:
    int j = 1;
    AnalysisDataLongFormatPlotModulePointer lfpltResMapping(new AnalysisDataLongFormatPlotModule(j));
    const char *fnResMapping = resMappingOutFileName_.c_str();
    std::vector<std::string> headerResMapping = {"t", 
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
    radiusFilePath.replace(radiusFilePath.begin() + lastSlash - 5, 
                           radiusFilePath.end(), 
                           "share/data/vdwradii/");
        
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
//    const Selection &initProbePosSelection = pdata -> parallelSelection(initProbePosSelection_);

    // get data handles for this frame:
    AnalysisDataHandle dhResMapping = pdata -> dataHandle(dataResMapping_);
    AnalysisDataHandle dhFrameStream = pdata -> dataHandle(frameStreamData_);

	// get data for frame number frnr into data handle:
    dhResMapping.startFrame(frnr, fr.time);
    dhFrameStream.startFrame(frnr, fr.time);


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
        pfm.reset(new InplaneOptimisedProbePathFinder(pfPar_,
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
        pfm.reset(new NaiveCylindricalPathFinder(pfPar_,
                                                 initProbePos,
                                                 chanDirVec));
    }

    // set parameters:
    pfm -> setParameters(pfParams_);




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
    std::cout.flush();
    clock_t tMolPath = std::clock();
    MolecularPath molPath = pfm -> getMolecularPath();
    tMolPath = (std::clock() - tMolPath)/CLOCKS_PER_SEC;
    std::cout<<"done in  "<<tMolPath<<" sec"<<std::endl;
    

    // which method do we use for path alignment?
    if( pfPathAlignmentMethod_ == ePathAlignmentMethodNone )
    {
        // no need to do anything in this case
    }
    else if( pfPathAlignmentMethod_ == ePathAlignmentMethodIpp )
    {
        // map initial probe position onto pathway:
        std::vector<gmx::RVec> ipp;
        ipp.push_back(initProbePos);
        std::vector<gmx::RVec> mappedIpp = molPath.mapPositions(
                ipp, 
                mappingParams_);

        // shift coordinates of molecular path appropriately:
        molPath.shift(mappedIpp.front());
    }

    // get original path points and radii:
    std::vector<gmx::RVec> pathPoints = molPath.pathPoints();
    std::vector<real> pathRadii = molPath.pathRadii();

    // add original path points to frame stream dataset:
    dhFrameStream.selectDataSet(1);
    for(size_t i = 0; i < pathPoints.size(); i++)
    {
        dhFrameStream.setPoint(0, pathPoints.at(i)[XX]);
        dhFrameStream.setPoint(1, pathPoints.at(i)[YY]);
        dhFrameStream.setPoint(2, pathPoints.at(i)[ZZ]);
        dhFrameStream.setPoint(3, pathRadii.at(i));
        dhFrameStream.finishPointSet();
    }

    // add radius spline knots and control points to frame stream dataset:
    dhFrameStream.selectDataSet(2);
    std::vector<real> radiusKnots = molPath.poreRadiusUniqueKnots();    
    std::vector<real> radiusCtrlPoints = molPath.poreRadiusCtrlPoints();
    for(size_t i = 0; i < radiusKnots.size(); i++)
    {
        dhFrameStream.setPoint(0, radiusKnots.at(i));
        dhFrameStream.setPoint(1, radiusCtrlPoints.at(i));
        dhFrameStream.finishPointSet();
    }
    
    // add centre line spline knots and control points to frame stream dataset:
    dhFrameStream.selectDataSet(3);
    std::vector<real> centreLineKnots = molPath.centreLineUniqueKnots();    
    std::vector<gmx::RVec> centreLineCtrlPoints = molPath.centreLineCtrlPoints();
    for(size_t i = 0; i < centreLineKnots.size(); i++)
    {
        dhFrameStream.setPoint(0, centreLineKnots.at(i));
        dhFrameStream.setPoint(1, centreLineCtrlPoints.at(i)[XX]);
        dhFrameStream.setPoint(2, centreLineCtrlPoints.at(i)[YY]);
        dhFrameStream.setPoint(3, centreLineCtrlPoints.at(i)[ZZ]);
        dhFrameStream.finishPointSet();
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
            mappingParams_);
    tMapResCog = (std::clock() - tMapResCog)/CLOCKS_PER_SEC;
    std::cout<<"mapped "<<poreCogMappedCoords.size()
             <<" particles in "<<1000*tMapResCog<<" ms"<<std::endl;

    // map pore residue C-alpha onto pathway:
    std::cout<<"mapping pore residue C-alpha onto pathway ... ";
    clock_t tMapResCal = std::clock();
    std::map<int, gmx::RVec> poreCalMappedCoords = molPath.mapSelection(
            poreMappingSelCal, 
            mappingParams_);
    tMapResCal = (std::clock() - tMapResCal)/CLOCKS_PER_SEC;
    std::cout<<"mapped "<<poreCalMappedCoords.size()
             <<" particles in "<<1000*tMapResCal<<" ms"<<std::endl;

    
    // check if particles are pore-lining:
    std::cout<<"checking which residues are pore-lining ... ";
    clock_t tResPoreLining = std::clock();
    std::map<int, bool> poreLining = molPath.checkIfInside(
            poreCogMappedCoords, 
            poreMappingMargin_);
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
    // TODO: this functionality should probably be handled outside the main analysis loop
    for(auto it = poreCogMappedCoords.begin(); it != poreCogMappedCoords.end(); it++)
    {
        SelectionPosition pos = poreMappingSelCog.position(it->first);
        
        // add points to dataset:
        dhResMapping.setPoint(0, pos.mappedId());                    // refId
        dhResMapping.setPoint(1, it -> second[SS]); // s
        dhResMapping.setPoint(2, it -> second[RR]); // rho
        dhResMapping.setPoint(3, it -> second[PP]);// phi
        dhResMapping.setPoint(4, poreLining[it -> first]);             // poreLining
        dhResMapping.setPoint(5, poreFacing[it -> first]);             // poreFacing TODO
        dhResMapping.finishPointSet();
    }
    
    // now add mapped residue coordinates to data handle:
    // FIXME JSON error caused here? --> only with cylindrical path finder!
    // --> nope, also with the other one if all legacy code if properly removed!
    // --> commenting this out certainly helps
    
    dhFrameStream.selectDataSet(4);
    
    // add mapped residues to data container:
    for(auto it = poreCogMappedCoords.begin(); it != poreCogMappedCoords.end(); it++)
    {
        dhFrameStream.setPoint(0, poreMappingSelCog.position(it -> first).mappedId());
        dhFrameStream.setPoint(1, it -> second[SS]);     // s
        dhFrameStream.setPoint(2, it -> second[RR]);     // rho
        dhFrameStream.setPoint(3, it -> second[PP]);     // phi
        dhFrameStream.setPoint(4, poreLining[it -> first]);     // pore lining?
        dhFrameStream.setPoint(5, poreFacing[it -> first]);     // pore facing?
        dhFrameStream.setPoint(6, poreMappingSelCog.position(it -> first).x()[0]);  // x
        dhFrameStream.setPoint(7, poreMappingSelCog.position(it -> first).x()[1]);  // y
        dhFrameStream.setPoint(8, poreMappingSelCog.position(it -> first).x()[2]);  // z
        dhFrameStream.finishPointSet();
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
            mappingParams_);
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
    dhFrameStream.selectDataSet(5);
    
    // add mapped residues to data container:
    for(auto it = solventMappedCoords.begin(); 
        it != solventMappedCoords.end(); 
        it++)
    {
         dhFrameStream.setPoint(0, solvMapSel.position(it -> first).mappedId()); // res.id
         dhFrameStream.setPoint(1, it -> second[0]);     // s
         dhFrameStream.setPoint(2, it -> second[1]);     // rho
         dhFrameStream.setPoint(3, 0.0);     // phi // FIXME wrong, but JSON cant handle NaN
         dhFrameStream.setPoint(4, solvInsidePore[it -> first]);     // inside pore
         dhFrameStream.setPoint(5, solvInsideSample[it -> first]);     // inside sample
         dhFrameStream.setPoint(6, solvMapSel.position(it -> first).x()[0]);  // x
         dhFrameStream.setPoint(7, solvMapSel.position(it -> first).x()[1]);  // y
         dhFrameStream.setPoint(8, solvMapSel.position(it -> first).x()[2]);  // z
         dhFrameStream.finishPointSet();
    }

    
    // ESTIMATE SOLVENT DENSITY
    //-------------------------------------------------------------------------

    // TODO this entire section can easily be made its own class

    // build a vector of sample points inside the pathway:
    std::vector<real> solventSampleCoordS;
    solventSampleCoordS.reserve(solventMappedCoords.size());
    for(auto isInsidePath : solvInsideSample)
    {
        // is this particle inside the pathway?
        if( isInsidePath.second )
        {
            // add arc length coordinate to sample vector:
            solventSampleCoordS.push_back(
                    solventMappedCoords[isInsidePath.first][SS]);
        }
    }

    // create density estimator:
    std::unique_ptr<AbstractDensityEstimator> densityEstimator;
    if( deMethod_ == eDensityEstimatorHistogram )
    {
        densityEstimator.reset(new HistogramDensityEstimator());
    }
    else if( deMethod_ == eDensityEstimatorKernel )
    {
        densityEstimator.reset(new KernelDensityEstimator());
    }

    // set parameters for density estimation:
    densityEstimator -> setParameters(deParams_);

    // estimate density of solvent particles along arc length coordinate:
    std::cout<<"estimating solvent density...";
    std::cout.flush();
    SplineCurve1D solventDensityCoordS = densityEstimator -> estimate(
            solventSampleCoordS);
    std::cout<<" done"<<std::endl;

    // add spline curve parameters to data handle:   
    dhFrameStream.selectDataSet(6);
    for(size_t i = 0; i < solventDensityCoordS.ctrlPoints().size(); i++)
    {
        dhFrameStream.setPoint(
                0, 
                solventDensityCoordS.uniqueKnots().at(i));
        dhFrameStream.setPoint(
                1, 
                solventDensityCoordS.ctrlPoints().at(i));
        dhFrameStream.finishPointSet();
    }

    // track range covered by solvent:
    real solventRangeLo = solventDensityCoordS.uniqueKnots().front();
    real solventRangeHi = solventDensityCoordS.uniqueKnots().back();


    // ADD AGGREGATE DATA TO PARALLELISABLE CONTAINER
    //-------------------------------------------------------------------------   

    // add aggegate path data:
    dhFrameStream.selectDataSet(0);

    // only one point per frame:
    dhFrameStream.setPoint(0, molPath.minRadius().second);
    dhFrameStream.setPoint(1, molPath.length());
    dhFrameStream.setPoint(2, molPath.volume());
    dhFrameStream.setPoint(3, numSolvInsidePore); 
    dhFrameStream.setPoint(4, numSolvInsideSample); 
    dhFrameStream.setPoint(5, solventRangeLo); 
    dhFrameStream.setPoint(6, solventRangeHi); 
    dhFrameStream.finishPointSet();


//    std::cout<<"solventRangeLo = "<<solventRangeLo<<"  ";
  //  std::cout<<"solventRangeHi = "<<solventRangeLo<<std::endl;


    // WRITE PORE TO OBJ FILE
    //-------------------------------------------------------------------------

    // TODO: this should be moved to a separate binary!

    MolecularPathObjExporter molPathExp;
    molPathExp(objOutputFileName_.c_str(),
               molPath);


    // FINISH FRAME
    //-------------------------------------------------------------------------

    std::cout<<std::endl;

	// finish analysis of current frame:
    dhResMapping.finishFrame();
    dhFrameStream.finishFrame();
}



/*
 *
 */
void
trajectoryAnalysis::finishAnalysis(int numFrames)
{
    // transfer file names from user input:
    std::string inFileName = std::string("stream_") + jsonOutputFileName_;
    std::string outFileName = jsonOutputFileName_;
    std::fstream inFile;
    std::fstream outFile;

    // READ PER-FRAME DATA AND AGGREGATE ALL NON-PROFILE DATA
    // ------------------------------------------------------------------------

    // TODO: this needs another pass through the infile and should collect
    // pore summary data like length and volume

    // openen per-frame data set for reading:
    inFile.open(inFileName, std::fstream::in);

    // prepare summary statistics for aggregate properties:
    SummaryStatistics minRadiusSummary;
    SummaryStatistics lengthSummary;
    SummaryStatistics volumeSummary;
    SummaryStatistics numPathSummary;
    SummaryStatistics numSampleSummary;
    SummaryStatistics solventRangeLoSummary;
    SummaryStatistics solventRangeHiSummary;

    // read file line by line and calculate summary statistics:
    int linesRead = 0;
    std::string line;
    while( std::getline(inFile, line) )
    {
        // read line into JSON document:
        rapidjson::StringStream lineStream(line.c_str());
        rapidjson::Document lineDoc;
        lineDoc.ParseStream(lineStream);

        // sanity checks:
        if( !lineDoc.IsObject() )
        {
            // FIXME this is where the JSON error occurs
            std::string error = "Line " + std::to_string(linesRead) + 
            " read from" + inFileName + " is not valid JSON object.";
            throw std::runtime_error(error);
        }
      
        // calculate summary statistics of aggregate variables:
        minRadiusSummary.update(
                lineDoc["pathSummary"]["minRadius"][0].GetDouble());
        lengthSummary.update(
                lineDoc["pathSummary"]["length"][0].GetDouble());
        volumeSummary.update(
                lineDoc["pathSummary"]["volume"][0].GetDouble());
        numPathSummary.update(
                lineDoc["pathSummary"]["numPath"][0].GetDouble());
        numSampleSummary.update(
                lineDoc["pathSummary"]["numSample"][0].GetDouble());
        solventRangeLoSummary.update(
                lineDoc["pathSummary"]["solventRangeLo"][0].GetDouble());
        solventRangeHiSummary.update(
                lineDoc["pathSummary"]["solventRangeHi"][0].GetDouble());

        // increment line counter:
        linesRead++;
    }

    // close per frame data set:
    inFile.close();
    
    // sanity check:
    if( linesRead != numFrames )
    {
        throw std::runtime_error("Number of frames read does not equal number"
        "of frames analyised.");
    }


    // READ PER-FRAME DATA AND AGGREGATE TIME-AVERAGED PORE PROFILE
    // ------------------------------------------------------------------------

    // define set of support points for profile evaluation:
    // FIXME this needs more than a heuristic!
    // FIXME also will not work when alignment = none is selected
    // TODO number of support points should be use settable
    std::vector<real> supportPoints;
    size_t numSupportPoints = 1000;
    real supportPointsLo = solventRangeLoSummary.min() + 2.0*deEvalRangeCutoff_*deBandWidth_;
    real supportPointsHi = solventRangeHiSummary.max() - 2.0*deEvalRangeCutoff_*deBandWidth_;
//    real supportPointsLo = -5.0;
//    real supportPointsHi = 5.0;
    real supportPointsStep = (supportPointsHi - supportPointsLo) / (numSupportPoints - 1);
    for(size_t i = 0; i < numSupportPoints; i++)
    {
        supportPoints.push_back(supportPointsLo + i*supportPointsStep);
    }
   
    /* // old solution:
    real extrapDist = mappingParams_.extrapDist_;
    real step = (lengthSummary.max() + 2.0*extrapDist) / (numSupportPoints - 1);
    for(size_t i = 0; i < numSupportPoints; i++)
    {
        supportPoints.push_back(-0.5*lengthSummary.max() - extrapDist + i*step);
    }
    */

    // open JSON data file in read mode:
    inFile.open(inFileName.c_str(), std::fstream::in);
    
    // prepare containers for profile summaries:
    std::vector<SummaryStatistics> radiusSummary(supportPoints.size());
    std::vector<SummaryStatistics> solventDensitySummary(supportPoints.size());
    std::vector<SummaryStatistics> energySummary(supportPoints.size());

    // read file line by line:
    int linesProcessed = 0;
    while( std::getline(inFile, line) )
    {
        std::cout<<"linesProcessed = "<<linesProcessed<<std::endl;

        // read line into JSON document:
        rapidjson::StringStream lineStream(line.c_str());
        rapidjson::Document lineDoc;
        lineDoc.ParseStream(lineStream);

        // sanity checks:
        if( !lineDoc.IsObject() )
        {
            std::string error = "Line " + std::to_string(linesProcessed) + 
            " read from" + inFileName + " is not valid JSON object.";
            throw std::runtime_error(error);
        }

        // create molecular path:
        MolecularPath molPath(lineDoc);

        // sample radius at support points and add to summary statistics:
        std::vector<real> radiusSample = molPath.sampleRadii(supportPoints); 
        for(size_t i = 0; i < radiusSample.size(); i++)
        {
            radiusSummary.at(i).update(radiusSample.at(i));
        }


        // TODO: this should get its own class:

        // get spline parameters from JSON:
        std::vector<real> solventDensityKnots;
        std::vector<real> solventDensityCtrlPoints;
        for(size_t i = 0; i < lineDoc["solventDensitySpline"]["knots"].Size(); i++)
        {
            solventDensityKnots.push_back(
                    lineDoc["solventDensitySpline"]["knots"][i].GetDouble());
            solventDensityCtrlPoints.push_back(
                    lineDoc["solventDensitySpline"]["ctrl"][i].GetDouble());
        }
        solventDensityKnots.push_back(
                solventDensityKnots.back());
        solventDensityKnots.insert(
                solventDensityKnots.begin(),
                solventDensityKnots.front());

        // construct Spline curve;
        SplineCurve1D solventDensitySpline(
                1,
                solventDensityKnots,
                solventDensityCtrlPoints);

        // sample from spline curve:
        std::vector<real> solventDensitySample;
        for(auto eval : supportPoints)
        {
            solventDensitySample.push_back(solventDensitySpline.evaluate(
                    eval,
                    0));
        }


        // get total number of particles in sample for this time step:
        int totalNumber = lineDoc["pathSummary"]["numSample"][0].GetDouble();

        // convert to number density and add to summary statistic:
        NumberDensityCalculator ndc;
        solventDensitySample = ndc(
                solventDensitySample, 
                radiusSample, 
                totalNumber);
        for(size_t i = 0; i < solventDensitySample.size(); i++)
        {
            solventDensitySummary.at(i).update(solventDensitySample.at(i));
        }
       
        // convert to energy and add to summary statistic:
        BoltzmannEnergyCalculator bec;
        std::vector<real> energySample = bec.calculate(solventDensitySample);
        for(size_t i = 0; i < energySample.size(); i++)
        {
            energySummary.at(i).update(energySample.at(i));
        }


        // increment line counter:
        linesProcessed++;
    }

    // sanity check:
    if( linesProcessed != numFrames )
    {
        std::string error = "Number of lines read from JSON file does not"
        "equal number of frames processed!"; 
        throw std::runtime_error(error);
    }

    // close filestream object:
    inFile.close();

    

    // CREATING OUTPUT JSON
    // ------------------------------------------------------------------------

    // TODO also need to write density and energy profiles

    // prepare JSON document for output:
    rapidjson::Document outDoc;
    rapidjson::Document::AllocatorType &alloc = outDoc.GetAllocator();
    outDoc.SetObject();

    // create JSON object for reproducibility information:
    // TODO this should probably get its own class
    rapidjson::Value reproInfo;
    reproInfo.SetObject();
    reproInfo.AddMember(
            "version",
            chapVersionString(),
            alloc);
    reproInfo.AddMember(
            "commandLine",
            std::string(gmx::getProgramContext().commandLine()),
            alloc);

    // add reproducibility information to output JSON:
    outDoc.AddMember(
            "reproducibilityInformation",
            reproInfo,
            alloc);

    // create JSON object for pore summary data:
    rapidjson::Value pathSummary;
    pathSummary.SetObject();

    SummaryStatisticsJsonConverter ssjc;
    pathSummary.AddMember(
            "minRadius",
            ssjc.convert(minRadiusSummary, outDoc.GetAllocator()),
            outDoc.GetAllocator());
    pathSummary.AddMember(
            "length",
            ssjc.convert(lengthSummary, outDoc.GetAllocator()),
            outDoc.GetAllocator());
    pathSummary.AddMember(
            "volume",
            ssjc.convert(volumeSummary, outDoc.GetAllocator()),
            outDoc.GetAllocator());
    pathSummary.AddMember(
            "numPath",
            ssjc.convert(numPathSummary, outDoc.GetAllocator()),
            outDoc.GetAllocator());
    pathSummary.AddMember(
            "numSample",
            ssjc.convert(numSampleSummary, outDoc.GetAllocator()),
            outDoc.GetAllocator());

    // add summary data to output document:
    outDoc.AddMember("pathSummary", pathSummary, alloc);


    // create JSON object for pore profile:
    rapidjson::Value pathProfile;
    pathProfile.SetObject();

    // create JSON arrays to hold pore profile values:
    rapidjson::Value supportPts(rapidjson::kArrayType);

    rapidjson::Value radiusMin(rapidjson::kArrayType);
    rapidjson::Value radiusMax(rapidjson::kArrayType);
    rapidjson::Value radiusMean(rapidjson::kArrayType);
    rapidjson::Value radiusSd(rapidjson::kArrayType);    

    rapidjson::Value densityMin(rapidjson::kArrayType);
    rapidjson::Value densityMax(rapidjson::kArrayType);
    rapidjson::Value densityMean(rapidjson::kArrayType);
    rapidjson::Value densitySd(rapidjson::kArrayType);    

    rapidjson::Value energyMin(rapidjson::kArrayType);
    rapidjson::Value energyMax(rapidjson::kArrayType);
    rapidjson::Value energyMean(rapidjson::kArrayType);
    rapidjson::Value energySd(rapidjson::kArrayType);    

    // fill JSON arrays with values::
    for(size_t i = 0; i < supportPoints.size(); i++)
    {
        // support points:
        supportPts.PushBack(supportPoints.at(i), alloc);

        // radius:
        radiusMin.PushBack(radiusSummary.at(i).min(), alloc);
        radiusMax.PushBack(radiusSummary.at(i).max(), alloc);
        radiusMean.PushBack(radiusSummary.at(i).mean(), alloc);
        radiusSd.PushBack(radiusSummary.at(i).sd(), alloc);

        // density:
        densityMin.PushBack(solventDensitySummary.at(i).min(), alloc);
        densityMax.PushBack(solventDensitySummary.at(i).max(), alloc);
        densityMean.PushBack(solventDensitySummary.at(i).mean(), alloc);
        densitySd.PushBack(solventDensitySummary.at(i).sd(), alloc);

        // energy:
        energyMin.PushBack(energySummary.at(i).min(), alloc);
        energyMax.PushBack(energySummary.at(i).max(), alloc);
        energyMean.PushBack(energySummary.at(i).mean(), alloc);
        energySd.PushBack(energySummary.at(i).sd(), alloc);
    }

    // add JSON arrays to pore profile object:
    pathProfile.AddMember("s", supportPts, alloc);

    pathProfile.AddMember("radiusMin", radiusMin, alloc);
    pathProfile.AddMember("radiusMax", radiusMax, alloc);
    pathProfile.AddMember("radiusMean", radiusMean, alloc);
    pathProfile.AddMember("radiusSd", radiusSd, alloc);

    pathProfile.AddMember("densityMin", densityMin, alloc);
    pathProfile.AddMember("densityMax", densityMax, alloc);
    pathProfile.AddMember("densityMean", densityMean, alloc);
    pathProfile.AddMember("densitySd", densitySd, alloc);

    pathProfile.AddMember("energyMin", energyMin, alloc);
    pathProfile.AddMember("energyMax", energyMax, alloc);
    pathProfile.AddMember("energyMean", energyMean, alloc);
    pathProfile.AddMember("energySd", energySd, alloc);
    
    // add pore profile to output document:
    outDoc.AddMember("pathProfile", pathProfile, alloc);


    // WRITING OUTPUT JSON TO FILE
    // ------------------------------------------------------------------------

    // stringify output document:
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    outDoc.Accept(writer);
    std::string outLine(buffer.GetString(), buffer.GetSize());
    
    // open outgoing file stream:
    outFile.open(outFileName, std::fstream::out);

    // write JSON output to file:
    outFile<<outLine<<std::endl;

    // close out file stream:
    outFile.close();


    // COPYING PER-FRAME DATA TO FINAL OUTPUT FILE
    // ------------------------------------------------------------------------
    
    // open file with per-frame data and output data:
    inFile.open(inFileName, std::fstream::in);
    outFile.open(outFileName, std::fstream::app);

    // append input file to output file line by line:    
    int linesCopied = 0;
    std::string copyLine;
    while( std::getline(inFile, copyLine) )
    {
        // append line to out file:
        outFile<<copyLine<<std::endl;

        // increment line counter:
        linesCopied++;
    }

    // close file streams:
    inFile.close();
    outFile.close();

    // sanity checks:
    if( linesCopied != numFrames )
    {
        throw std::runtime_error("Could not copy all lines from per-frame data"
        "file to output data file.");
    }

    // delete temporary file:
    std::remove(inFileName.c_str());
}




void
trajectoryAnalysis::writeOutput()
{

}

