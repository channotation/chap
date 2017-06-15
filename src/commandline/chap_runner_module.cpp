#include <gromacs/options/timeunitmanager.h>
#include <gromacs/pbcutil/pbc.h>

#include "commandline/chap_runner_module.hpp"
#include "commandline/paralleloptions.h"


/*
 *
 */
void
ChapRunnerModule::init(gmx::CommandLineModuleSettings *settings)
{

}


/*
 *
 */
void
ChapRunnerModule::initOptions(gmx::IOptionsContainer *options,
                              gmx::ICommandLineOptionsModuleSettings *settings)
{
    std::cerr<<"BEGIN INIT OPTIONS"<<std::endl;

    std::shared_ptr<gmx::TimeUnitBehavior> timeUnitBehavior(
        new gmx::TimeUnitBehavior());
    std::shared_ptr<gmx::SelectionOptionBehavior> selectionOptionBehavior(
        new gmx::SelectionOptionBehavior(&selections_, &topologyProvider_));
    settings -> addOptionsBehavior(timeUnitBehavior);
    settings -> addOptionsBehavior(selectionOptionBehavior);
   
    gmx::IOptionsContainer &commonOptions = options -> addGroup();
    gmx::IOptionsContainer &moduleOptions = options -> addGroup();

    settings_.setOptionsModuleSettings(settings);
    module_ -> initOptions(&moduleOptions, &settings_);
    common_.initOptions(&commonOptions, timeUnitBehavior.get());
    settings_.setOptionsModuleSettings(nullptr);
    selectionOptionBehavior -> initOptions(&commonOptions);

    std::cerr<<"END INIT OPTIONS"<<std::endl;
}


/*
 *
 */
void
ChapRunnerModule::optionsFinished()
{
    std::cerr<<"BEGIN OPTIONS FINISHED"<<std::endl;

    common_.optionsFinished();
    module_ -> optionsFinished(&settings_);

    std::cerr<<"END OPTIONS FINISHED"<<std::endl;
}







/*
 *
 */
int
ChapRunnerModule::run()
{
    std::cerr<<"BEGIN RUN"<<std::endl;
 

    // initialise analysis:
    //-------------------------------------------------------------------------

    common_.initTopology();
    const ChapTopologyInformation &topology = common_.topologyInformation();
    module_ -> initAnalysis(settings_, topology);


    // load first frame:
    //-------------------------------------------------------------------------

    common_.initFirstFrame();
    common_.initFrameIndexGroup();
    module_ -> initAfterFirstFrame(settings_, common_.frame());

    t_pbc pbc;
    t_pbc *ppbc = settings_.hasPBC() ? &pbc : NULL;


    // loop over frames:
    //-------------------------------------------------------------------------


std::cout<<"loop over frames"<<std::endl;

    int nFrames = 0;

    gmx::AnalysisDataParallelOptions dataOptions;
    gmx::TrajectoryAnalysisModuleDataPointer pdata(
            module_ -> startFrames(dataOptions, selections_));


    // actual loop over frames:
    // TODO: parallelise this
    do
    {
        // read and access trajectory frame:
        common_.initFrame();
        t_trxframe &frame = common_.frame();

        // set periodic boundary conditions:
        if( ppbc != NULL )
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        // reevaluate all selections for this frame:
        selections_.evaluate(&frame, ppbc);

        // perform analysis on this frame:
        module_ -> analyzeFrame(nFrames, frame, ppbc, pdata.get());
        module_ -> finishFrameSerial(nFrames);

        // increment frame counter:
        ++nFrames;
    }
    while( common_.readNextFrame() );
    module_ -> finishFrames(pdata.get());


    // finish data:
    if( pdata.get() != NULL )
    {
        pdata->finish();
    }
    pdata.reset();

    // inform user:
    if( common_.hasTrajectory() )
    {
        fprintf(stdout, 
                "Analysed %d frames, last time %.f\n",
                nFrames,
                common_.frame().time);
    }
    else
    {
        fprintf(stdout, "Analysed topology coordinates\n");
    }

    // restore dynamic selections:
    selections_.evaluateFinal(nFrames);

    // finish the analysis:
    module_ -> finishAnalysis(nFrames);
    module_ -> writeOutput();


    std::cerr<<"END RUN"<<std::endl;

    return 0;

}




















