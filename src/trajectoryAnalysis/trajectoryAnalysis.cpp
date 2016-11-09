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

	   /* 
    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnDist_).defaultBasename("avedist")
                           .description("Average distances from reference group"));
                           .description("Cutoff for distance calculation (0 = no cutoff)"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
	*/
}




/*
 * 
 */
void
trajectoryAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                                 const TopologyInformation         & /*top*/)
{
	// set cutoff distance for grid search as specified in user input:
	nb_.setCutoff(cutoff_);
	std::cout<<"Setting cutoff to: "<<cutoff_<<std::endl;

	// set number of columns in data set (one column per small particle type):
	data_.setColumnCount(0, sel_.size());
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

	// initialise neighbourhood search:
	AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, ref_selection);


	// loop over small particle selections:
	for(size_t g = 0; g < sel_.size(); ++g)
	{
		// get thread-local selection of small particles:
		const Selection &sel = pdata -> parallelSelection(sel_[g]);

		// get number of particles in selection:
		int n_part = sel.posCount();
		
	}







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





