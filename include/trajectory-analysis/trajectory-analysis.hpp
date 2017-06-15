#ifndef TRAJECTORYANALYSIS_HPP
#define TRAJECTORYANALYSIS_HPP

#include <iostream>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/vdw_radius_provider.hpp"

#include "commandline/chap_traj_ana_runner_common.hpp"

using namespace gmx;


//
class ChapTopologyInformation;


/*
 *
 */
class ChapTrajectoryAnalysisModule : public gmx::TrajectoryAnalysisModule
{
    public:

        ~ChapTrajectoryAnalysisModule(){};

        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const ChapTopologyInformation &top) = 0;

};

/*
 *
 */
typedef std::unique_ptr<ChapTrajectoryAnalysisModule> ChapTrajectoryAnalysisModulePointer;


class trajectoryAnalysis : public ChapTrajectoryAnalysisModule
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

    //
    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const ChapTopologyInformation &top){};
	
	// ??
	virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);
	
	// ??
	virtual void finishAnalysis(int nframes);
	
	// ??
	virtual void writeOutput();


    private:

    std::string poreParticleFileName_;  // positions of probe particles
    std::string smallParticleFileName_; // positions of small particles (e.g. water)
    std::string poreProfileFileName_;       // time averaged radius, energy etc.
    bool poreFile_;

	class ModuleData;
	double                           cutoff_;		// cutoff for grid search
	Selection                        refsel_;   	// selection of the reference group
	Selection                        ippsel_;   	// selection of the initial probe position group
    bool                             ippselIsSet_;
	SelectionList                    sel_;			// selection of the small particle groups
	AnalysisNeighborhood             nb_;			// neighbourhood for grid searches

    AnalysisData                     data_;			// raw data container
    AnalysisData                     dataResMapping_;


	std::unordered_map<int, real>	 vdwRadii_;		// vdW radii of all atoms
	real 							 maxVdwRadius_;	// largest vdW radius of all atoms


    int                              nOutPoints_;   // number of points on path sample



    // selection and topology for initial probe position:
    gmx::SelectionCollection initProbePosCollection_;
    gmx::Selection initProbePosSelection_;

    // path finding method parameters:
    real pfDefaultVdwRadius_;
    bool pfDefaultVdwRadiusIsSet_;
    eVdwRadiusDatabase pfVdwRadiusDatabase_;
    std::string pfVdwRadiusJson_;
    bool pfVdwRadiusJsonIsSet_;
    std::string pfMethod_;
    real pfProbeStepLength_;
    real pfProbeRadius_;
    real pfMaxFreeDist_;
    int pfMaxProbeSteps_;
    std::vector<real> pfInitProbePos_;
    bool pfInitProbePosIsSet_;
    std::vector<real> pfChanDirVec_;
    bool pfChanDirVecIsSet_;

    // simulated annealing parameters:
    int saRandomSeed_;
    int saMaxCoolingIter_;
    int saNumCostSamples_;
    real saXi_;
    real saConvRelTol_;
    real saInitTemp_;
    real saCoolingFactor_;
    real saStepLengthFactor_;
    bool saUseAdaptiveCandGen_;

    // Nelder-Mead parameters:
    int nmMaxIter_;





    bool debug_output_;

    // map for path finding parameters:
    std::map<std::string, real> pfParams_;

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

