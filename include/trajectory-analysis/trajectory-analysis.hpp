#ifndef TRAJECTORYANALYSIS_HPP
#define TRAJECTORYANALYSIS_HPP

#include <iostream>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

using namespace gmx;




class trajectoryAnalysis : public TrajectoryAnalysisModule
{
    public:

	// constructor for the trajectoryAnalysis module:
	trajectoryAnalysis();

	// method for adding all options of the trajectoryAnalysis module:
	virtual void initOptions(IOptionsContainer *options,
							 TrajectoryAnalysisSettings *settings);
	
	// ??
	virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation &top);
	
	// ??
	virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);
	
	// ??
	virtual void finishAnalysis(int nframes);
	
	// ??
	virtual void writeOutput();


    private:

	class ModuleData;
	double                           cutoff_;		// cutoff for grid search
	Selection                        refsel_;   	// selection of the reference group
	Selection                        ippsel_;   	// selection of the initial probe position group
    bool                             ippselIsSet_;
	SelectionList                    sel_;			// selection of the small particle groups
	AnalysisNeighborhood             nb_;			// neighbourhood for grid searches
	AnalysisData                     data_;			// raw data container
	std::vector<real>				 vdwRadii_;		// vdW radii of all atoms
	real 							 maxVdwRadius_;	// largest vdW radius of all atoms

    std::string pfMethod_;

    real pfProbeStepLength_;
    real pfProbeRadius_;
    real pfMaxFreeDist_;

    int pfMaxProbeSteps_;

    std::vector<real> pfInitProbePos_;
    bool pfInitProbePosIsSet_;
    std::vector<real> pfChanDirVec_;
    bool pfChanDirVecIsSet_;

    int saRandomSeed_;
    int saMaxCoolingIter_;
    int saNumCostSamples_;

    real saXi_;
    real saConvRelTol_;
    real saInitTemp_;
    real saCoolingFactor_;
    real saStepLengthFactor_;

    bool saUseAdaptiveCandGen_;

    bool debug_output_;

	// calculate the radius of a spherical void with given centre: 
	real calculateVoidRadius(RVec centre,
                             t_pbc *pbc,
							 const Selection refSelection);

	// optimise centre coordinates for maximum void radius:
	real maximiseVoidRadius(RVec &centre,
							RVec channelVec,
							t_pbc *pbc,
							const Selection refSelection);
};


#endif

