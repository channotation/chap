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
	SelectionList                    sel_;			// selection of the small particle groups
	AnalysisNeighborhood             nb_;			// neighbourhood for grid searches
	AnalysisData                     data_;			// raw data container
	std::vector<real>				 vdwRadii;		// vdW radii of all atoms
	real 							 maxVdwRadius;	// largest vdW radius of all atoms

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

