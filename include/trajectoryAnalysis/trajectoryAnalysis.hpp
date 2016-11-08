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
        trajectoryAnalysis();
        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;
        std::string                      fnDist_;
        double                           cutoff_;
        Selection                        refsel_;
        SelectionList                    sel_;
        AnalysisNeighborhood             nb_;
        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};


#endif

