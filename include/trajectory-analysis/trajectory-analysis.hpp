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

#include "io/pdb_io.hpp"

#include "path-finding/abstract_path_finder.hpp"
#include "path-finding/molecular_path.hpp"
#include "path-finding/vdw_radius_provider.hpp"

#include "statistics/abstract_density_estimator.hpp"

using namespace gmx;



class ChapTrajectoryAnalysis : public TrajectoryAnalysisModule
{
    public:

        // constructor for the trajectoryAnalysis module:
        ChapTrajectoryAnalysis();

        // method for adding all options of the trajectoryAnalysis module:
        virtual void initOptions(
                IOptionsContainer *options,
                TrajectoryAnalysisSettings *settings);
        
        //
        virtual void initAnalysis(
                const TrajectoryAnalysisSettings &settings,
                const TopologyInformation &top);

        //
        virtual void initAfterFirstFrame(
                const TrajectoryAnalysisSettings &settings,
                const t_trxframe &fr);


        //
        virtual void analyzeFrame(
                int frnr, 
                const t_trxframe &fr, 
                t_pbc *pbc,
                TrajectoryAnalysisModuleData *pdata);
        
        // ??
        virtual void finishAnalysis(int nframes);
        
        // ??
        virtual void writeOutput();


    private:

        // find file path fo index files:
        virtual void obtainNdxFilePathInfo();   
        std::string customNdxFileName_;

        
        // check input parameter validity:
        virtual void checkParameters();

        
        // names of output files:
        std::string outputBaseFileName_;
        std::string outputJsonFileName_;
        std::string outputPdbFileName_;

        
        // selections:
        SelectionList sel_;	    // selection of the small particle groups
        Selection refsel_;   	// selection of the reference group
        Selection ippsel_;   	// selection of the initial probe position group
        bool ippselIsSet_;

        
        // internal selections for pore mapping:
        std::string pfSelString_;
        SelectionCollection poreMappingSelCol_;
        SelectionCollection solvMappingSelCol_;
        Selection poreMappingSelCal_;
        Selection poreMappingSelCog_;
        Selection solvMappingSelCog_;
        real poreMappingMargin_;
        bool findPfResidues_;


        // data containers:
        AnalysisData frameStreamData_;
        AnalysisData dataResMappingPdb_;
        AnalysisData timingData_; 


        // pore residue chemical and physical information:
        eHydrophobicityDatabase hydrophobicityDatabase_;
        bool hydrophobicityDatabaseIsSet_;
        real hydrophobicityDefault_;
        bool hydrophobicityDefaultIsSet_;
        std::string hydrophobicityJson_;
        bool hydrophobicityJsonIsSet_;
        ResidueInformationProvider resInfo_;
        

        // output parameters:
        int outputNumPoints_;        
        real outputExtrapDist_;
        real outputGridSampleDist_;
        real outputCorrectionThreshold_;
        bool outputDetailed_;
        PdbStructure outputStructure_;


        // path finding:
        double cutoff_;
        bool cutoffIsSet_;
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
        std::map<std::string, real> pfPar_;
        std::unordered_map<int, real> vdwRadii_;
        real maxVdwRadius_;


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

        
        // density estimation parameters:
        eDensityEstimator deMethod_;
        DensityEstimationParameters deParams_;
        real deResolution_;
        real deBandWidth_;
        real deBandWidthScale_;
        real deEvalRangeCutoff_;


        // hydrophobicity profile parameters:
        real hpBandWidth_;
        real hpEvalRangeCutoff_;
        real hpResolution_;
        DensityEstimationParameters hydrophobKernelParams_;
        
        
        // molecular pathway for first frame:
        std::unique_ptr<MolecularPath> molPathAvg_;
};

#endif

