#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h" 
#include "gromacs/utility/smalloc.h" 
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/programcontext.h"

#include "commandline/chap_traj_ana_runner_common.hpp"
#include "commandline/chap_topology_information.hpp"

#include <iostream>


/******************************************************************************
 * RUNNER COMMON IMPL
 *****************************************************************************/


/*
 *
 */
class ChapTrajAnaRunnerCommon::Impl : public gmx::ITopologyProvider
{
    public:

        explicit Impl(gmx::TrajectoryAnalysisSettings *settings);
        ~Impl();

        bool hasTrajectory() const { return true; };

        void initTopology(bool required);
        void initFirstFrame();
        void initFrameIndexGroup();
        void finishTrajectory();

        // implement interface of ITopologyProvider:
        virtual t_topology *getTopology(bool required)
        {
            initTopology(required);
            return topInfo_.topology();
        }
        virtual int getAtomCount()
        {
            if( !topInfo_.hasTopology() )
            {
                if( trajectoryGroup_.isValid() )
                {
                    GMX_THROW(gmx::InconsistentInputError("-fgroup is only supported if -s is also specified"));
                }
                initFirstFrame();
                return fr -> natoms;
            }
            return -1;
        }

        gmx::TrajectoryAnalysisSettings &settings_; 
        ChapTopologyInformation topInfo_; 

        std::string trjfile_; 
        std::string topfile_;  

        gmx::Selection trajectoryGroup_;  
        double startTime_; 
        double endTime_;  
        double deltaTime_;
        bool bStartTimeSet_; 
        bool bEndTimeSet_; 
        bool bDeltaTimeSet_;
        bool bTrajOpen_;
        t_trxframe *fr; 
        gmx_rmpbc_t gpbc_; 
        t_trxstatus *status_; 
        gmx_output_env_t *oenv_; 
};


/*
 *
 */
ChapTrajAnaRunnerCommon::Impl::Impl(gmx::TrajectoryAnalysisSettings *settings)
    : settings_(*settings)
    , startTime_(0.0)
    , endTime_(0.0)
    , deltaTime_(0.0)
    , bStartTimeSet_(false)
    , bEndTimeSet_(false)
    , bDeltaTimeSet_(false)
    , bTrajOpen_(false)
    , fr(NULL)
    , gpbc_(NULL)
    , status_(NULL)
    , oenv_(NULL)
{

}



/*
 *
 */
ChapTrajAnaRunnerCommon::Impl::~Impl()
{
    finishTrajectory();
    if( fr != nullptr )
    {
        sfree(fr -> x);
        sfree(fr -> v);
        sfree(fr -> f);
        sfree(fr -> index);
        sfree(fr);
    }
    if( oenv_ != nullptr )
    {
        output_env_done(oenv_);
    }
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::Impl::initTopology(bool required)
{
    std::cout<<"COMMON: BEGIN INIT TOPOLOGY"<<std::endl;

    // Return immediately if the topology has already been loaded.
    if (topInfo_.hasTopology())
    {
        return;
    }

    if (required && topfile_.empty())
    {
        GMX_THROW(gmx::InconsistentInputError("No topology provided, but one is required for analysis"));
    }

    // Load the topology if requested.
    if (!topfile_.empty())
    {
        snew(topInfo_.top_, 1);
        topInfo_.bTop_ = read_tps_conf(topfile_.c_str(), topInfo_.top_, &topInfo_.ePBC_,
                                       &topInfo_.xtop_, NULL, topInfo_.boxtop_, TRUE);
        if (hasTrajectory()
            && !settings_.hasFlag(gmx::TrajectoryAnalysisSettings::efUseTopX))
        {
            sfree(topInfo_.xtop_);
            topInfo_.xtop_ = NULL;
        }
    }

    std::cout<<"COMMON: BEGIN INIT TOPOLOGY"<<std::endl;    
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::Impl::initFirstFrame()
{
    std::cout<<"COMMON: BEGIN INIT FIRST FRAME"<<std::endl;

    // Return if we have already initialized the trajectory.
    if (fr != NULL)
    {
        return;
    }
    time_unit_t time_unit
        = static_cast<time_unit_t>(settings_.timeUnit() + 1);
    output_env_init(&oenv_, gmx::getProgramContext(), time_unit, FALSE, exvgNONE, 0);

    int frflags = settings_.frflags();
    frflags |= TRX_NEED_X;

    snew(fr, 1);

    if (hasTrajectory())
    {
        if (!read_first_frame(oenv_, &status_, trjfile_.c_str(), fr, frflags))
        {
            GMX_THROW(gmx::FileIOError("Could not read coordinates from trajectory"));
        }
        bTrajOpen_ = true;

        if (topInfo_.hasTopology())
        {
            const int topologyAtomCount = topInfo_.topology()->atoms.nr;
            if (fr->natoms > topologyAtomCount)
            {
                const std::string message
                    = gmx::formatString("Trajectory (%d atoms) does not match topology (%d atoms)",
                                   fr->natoms, topologyAtomCount);
                GMX_THROW(gmx::InconsistentInputError(message));
            }
        }
    }
    else
    {
        // Prepare a frame from topology information.
        // TODO: Initialize more of the fields.
        if (frflags & (TRX_NEED_V))
        {
            GMX_THROW(gmx::NotImplementedError("Velocity reading from a topology not implemented"));
        }
        if (frflags & (TRX_NEED_F))
        {
            GMX_THROW(gmx::InvalidInputError("Forces cannot be read from a topology"));
        }
        fr->natoms = topInfo_.topology()->atoms.nr;
        fr->bX     = TRUE;
        snew(fr->x, fr->natoms);
        memcpy(fr->x, topInfo_.xtop_,
               sizeof(*fr->x) * fr->natoms);
        fr->bBox   = TRUE;
        copy_mat(topInfo_.boxtop_, fr->box);
    }

    set_trxframe_ePBC(fr, topInfo_.ePBC());
    if (topInfo_.hasTopology() && settings_.hasRmPBC())
    {
        gpbc_ = gmx_rmpbc_init(&topInfo_.topology()->idef, topInfo_.ePBC(),
                               fr->natoms);
    }
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::Impl::initFrameIndexGroup()
{
    if (!trajectoryGroup_.isValid())
    {
        return;
    }
    GMX_RELEASE_ASSERT(bTrajOpen_,
                       "Trajectory index only makes sense with a real trajectory");
    if (trajectoryGroup_.atomCount() != fr->natoms)
    {
        const std::string message = gmx::formatString(
                    "Selection specified with -fgroup has %d atoms, but "
                    "the trajectory (-f) has %d atoms.",
                    trajectoryGroup_.atomCount(), fr->natoms);
        GMX_THROW(gmx::InconsistentInputError(message));
    }
    fr->bIndex = TRUE;
    snew(fr->index, trajectoryGroup_.atomCount());
    std::copy(trajectoryGroup_.atomIndices().begin(),
              trajectoryGroup_.atomIndices().end(),
              fr->index);

}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::Impl::finishTrajectory()
{
    if (bTrajOpen_)
    {
        close_trx(status_);
        bTrajOpen_ = false;
    }
    if (gpbc_ != NULL)
    {
        gmx_rmpbc_done(gpbc_);
        gpbc_ = NULL;
    }
}


/******************************************************************************
 * RUNNER COMMON
 *****************************************************************************/

/*
 *
 */
ChapTrajAnaRunnerCommon::ChapTrajAnaRunnerCommon(
        gmx::TrajectoryAnalysisSettings *settings)
    : impl_(new Impl(settings))
{

}


/*
 *
 */
ChapTrajAnaRunnerCommon::~ChapTrajAnaRunnerCommon()
{

}


/*
 *
 */
gmx::ITopologyProvider *
ChapTrajAnaRunnerCommon::topologyProvider()
{
    std::cerr<<"BEGIN topologyProvider"<<std::endl;

    return impl_.get();

    std::cerr<<"END topologyProvider"<<std::endl;
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::initOptions(gmx::IOptionsContainer *options,
                                     gmx::TimeUnitBehavior *timeUnitBehavior)
{
    gmx::TrajectoryAnalysisSettings &settings = impl_->settings_;

    // Add common file name arguments.
    options->addOption(gmx::FileNameOption("f")
                           .filetype(gmx::eftTrajectory).inputFile()
                           .store(&impl_->trjfile_)
                           .defaultBasename("traj")
                           .description("Input trajectory or single configuration"));
    options->addOption(gmx::FileNameOption("s")
                           .filetype(gmx::eftTopology).inputFile()
                           .store(&impl_->topfile_)
                           .defaultBasename("topol")
                           .description("Input structure"));

    // Add options for trajectory time control.
    options->addOption(gmx::DoubleOption("b")
                           .store(&impl_->startTime_).storeIsSet(&impl_->bStartTimeSet_)
                           .timeValue()
                           .description("First frame (%t) to read from trajectory"));
    options->addOption(gmx::DoubleOption("e")
                           .store(&impl_->endTime_).storeIsSet(&impl_->bEndTimeSet_)
                           .timeValue()
                           .description("Last frame (%t) to read from trajectory"));
    options->addOption(gmx::DoubleOption("dt")
                           .store(&impl_->deltaTime_).storeIsSet(&impl_->bDeltaTimeSet_)
                           .timeValue()
                           .description("Only use frame if t MOD dt == first time (%t)"));

/*
    // Add time unit option.
    timeUnitBehavior->setTimeUnitFromEnvironment();
    timeUnitBehavior->addTimeUnitOption(options, "tu");
    timeUnitBehavior->setTimeUnitStore(&impl_->settings_.impl_->timeUnit);

    options->addOption(SelectionOption("fgroup")
                           .store(&impl_->trajectoryGroup_)
                           .onlySortedAtoms().onlyStatic()
                           .description("Atoms stored in the trajectory file "
                                        "(if not set, assume first N atoms)"));

    // Add plot options.
    settings.impl_->plotSettings.initOptions(options);

    // Add common options for trajectory processing.
    if(!settings.hasFlag( gmx::TrajectoryAnalysisSettings::efNoUserRmPBC) )
    {
        options->addOption(BooleanOption("rmpbc").store(&settings.impl_->bRmPBC)
                               .description("Make molecules whole for each frame"));
    }
    if( !settings.hasFlag(gmx::TrajectoryAnalysisSettings::efNoUserPBC) )
    {
        options->addOption(BooleanOption("pbc").store(&settings.impl_->bPBC)
                               .description("Use periodic boundary conditions for distance calculation"));
    }
    */
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::optionsFinished()
{
    if (impl_->trjfile_.empty() && impl_->topfile_.empty())
    {
        GMX_THROW(gmx::InconsistentInputError("No trajectory or topology provided, nothing to do!"));
    }

    if (impl_->trajectoryGroup_.isValid() && impl_->trjfile_.empty())
    {
        GMX_THROW(gmx::InconsistentInputError("-fgroup only makes sense together with a trajectory (-f)"));
    }
/*
    impl_->settings_.impl_->plotSettings.setTimeUnit(impl_->settings_.timeUnit());

    if (impl_->bStartTimeSet_)
    {
        setTimeValue(TBEGIN, impl_->startTime_);
    }
    if (impl_->bEndTimeSet_)
    {
        setTimeValue(TEND, impl_->endTime_);
    }
    if (impl_->bDeltaTimeSet_)
    {
        setTimeValue(TDELTA, impl_->deltaTime_);
    }    
    */
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::initTopology()
{
    const bool topologyRequired =
        impl_->settings_.hasFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    impl_->initTopology(topologyRequired);
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::initFirstFrame()
{
    impl_->initFirstFrame();
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::initFrameIndexGroup()
{
    impl_->initFrameIndexGroup();
}


/*
 *
 */
bool
ChapTrajAnaRunnerCommon::readNextFrame()
{
    bool bContinue = false;
    if (hasTrajectory())
    {
        bContinue = read_next_frame(impl_->oenv_, impl_->status_, impl_->fr);
    }
    if (!bContinue)
    {
        impl_->finishTrajectory();
    }
    return bContinue;
}


/*
 *
 */
void
ChapTrajAnaRunnerCommon::initFrame()
{
    if (impl_->gpbc_ != NULL)
    {
        gmx_rmpbc_trxfr(impl_->gpbc_, impl_->fr);
    }
}


/*
 *
 */
bool
ChapTrajAnaRunnerCommon::hasTrajectory() const
{
    return impl_ -> hasTrajectory();
}


/*
 *
 */
const ChapTopologyInformation &
ChapTrajAnaRunnerCommon::topologyInformation() const
{
    return impl_->topInfo_;
}


/*
 *
 */
t_trxframe &
ChapTrajAnaRunnerCommon::frame() const
{
    GMX_RELEASE_ASSERT(impl_->fr != NULL, "Frame not available when accessed");
    return *impl_->fr;
}

