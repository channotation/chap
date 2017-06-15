#ifndef CHAP_TRAJ_ANA_RUNNER_COMMON_HPP
#define CHAP_TRAJ_ANA_RUNNER_COMMON_HPP

#include <gromacs/selection/selectionoptionbehavior.h>
#include <gromacs/options/timeunitmanager.h> 
#include <gromacs/trajectoryanalysis.h>
#include "gromacs/fileio/trxio.h"
#include "gromacs/pbcutil/rmpbc.h"

#include "commandline/chap_topology_information.hpp"


/*
 *
 */
class ChapTrajAnaRunnerCommon
{
    public:

        // constructor and destructor:
        explicit ChapTrajAnaRunnerCommon(gmx::TrajectoryAnalysisSettings *settings);
        ~ChapTrajAnaRunnerCommon();

        //
        gmx::ITopologyProvider *topologyProvider();

        // 
        void initOptions(gmx::IOptionsContainer *options, 
                         gmx::TimeUnitBehavior *timeUnitBehavior);

        void optionsFinished();


        void initTopology();


        void initFirstFrame();


        void initFrameIndexGroup();


        bool readNextFrame();


        void initFrame();


        bool hasTrajectory() const;


        const ChapTopologyInformation &topologyInformation() const;


        t_trxframe &frame() const;


    private:

        class Impl;

        gmx::PrivateImplPointer<Impl> impl_;


};



#endif

