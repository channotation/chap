// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef TRAJECTORYANALYSIS_HPP
#define TRAJECTORYANALYSIS_HPP

#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_map>
#include <string>

#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "analysis-setup/residue_information_provider.hpp"

#include "path-finding/abstract_path_finder.hpp"
#include "path-finding/molecular_path.hpp"
#include "path-finding/vdw_radius_provider.hpp"

#include "statistics/abstract_density_estimator.hpp"

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
//	virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
//                              const TopologyInformation &top);

    //
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

    std::string jsonOutputFileName_;
    std::string objOutputFileName_;
    std::string resMappingOutFileName_;
    bool poreFile_;

	class ModuleData;
	double                           cutoff_;		// cutoff for grid search
    bool                             cutoffIsSet_;
	Selection                        refsel_;   	// selection of the reference group
	Selection                        ippsel_;   	// selection of the initial probe position group

    // internal selections for pore mapping:
    SelectionCollection              poreMappingSelCol_;
    SelectionCollection              solvMappingSelCol_;
    Selection                        poreMappingSelCal_;
    Selection                        poreMappingSelCog_;
    Selection                        solvMappingSelCog_;
    real                             poreMappingMargin_;



    bool                             ippselIsSet_;
	SelectionList                    sel_;			// selection of the small particle groups
	AnalysisNeighborhood             nb_;			// neighbourhood for grid searches

    AnalysisData                     frameStreamData_;
    AnalysisData                     dataResMapping_;
    AnalysisData                     dataResMappingPdb_;


	std::unordered_map<int, real>	 vdwRadii_;		// vdW radii of all atoms
	real 							 maxVdwRadius_;	// largest vdW radius of all atoms


    // pore residue chemical and physical information:
    eHydrophobicityDatabase hydrophobicityDatabase_;
    bool hydrophobicityDatabaseIsSet_;
    real hydrophobicityDefault_;
    bool hydrophobicityDefaultIsSet_;
    std::string hydrophobicityJson_;
    bool hydrophobicityJsonIsSet_;
    ResidueInformationProvider resInfo_;
    
    // hydrophobicity profile parameters:
    real hpBandWidth_;
    real hpEvalRangeCutoff_;
    real hpResolution_;
    DensityEstimationParameters hydrophobKernelParams_;


    // pore particle and group indices:
    std::vector<int> poreCAlphaIndices_;                    // c-alpha atomIds
    std::vector<int> residueIndices_;                       // all resIds
    std::vector<int> poreResidueIndices_;                   // resIds of pore selection
    std::vector<int> poreAtomIndices_;                      // atomIds of pore selection
    std::map<int, int> atomResidueMapping_;                 // maps atomId onto resId
    std::map<int, std::vector<int>> residueAtomMapping_;    // maps resId onto atomId



    int nOutPoints_;   // number of points on path sample



    // selection and topology for initial probe position:
    gmx::SelectionCollection initProbePosCollection_;
    gmx::Selection initProbePosSelection_;

    // path finding method parameters:
    real pfDefaultVdwRadius_;
    bool pfDefaultVdwRadiusIsSet_;
    eVdwRadiusDatabase pfVdwRadiusDatabase_;
    std::string pfVdwRadiusJson_;
    bool pfVdwRadiusJsonIsSet_;
    ePathFindingMethod pfMethod_;
    real pfProbeStepLength_;
    real pfProbeRadius_;
    real pfMaxProbeRadius_;
    int pfMaxProbeSteps_;
    std::vector<real> pfInitProbePos_;
    bool pfInitProbePosIsSet_;
    std::vector<real> pfChanDirVec_;
    bool pfChanDirVecIsSet_;
    ePathAlignmentMethod pfPathAlignmentMethod_;
    PathFindingParameters pfParams_;


    // simulated annealing parameters:
    int64_t saRandomSeed_;
    bool saRandomSeedIsSet_;
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

    // path mapping parameters:
    PathMappingParameters mappingParams_;

    // density estimation parameters:
    eDensityEstimator deMethod_;
    DensityEstimationParameters deParams_;
    real deResolution_;
    real deBandWidth_;
    real deEvalRangeCutoff_;


    bool debug_output_;

    // map for path finding parameters:
    std::map<std::string, real> pfPar_;

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

