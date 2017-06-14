#include "commandline/chap_runner_module.hpp"


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
/*
    // TODO: topology provider
    std::shared_ptr<gmx::TimeUnitBehaviour> timeUnitBehaviour(
        new gmx::TimeUnitBehaviour());
//   std::shared_ptr<SelectionOptionBehaviour> selectionOptionBehaviour(
//       new SelectioNoptionBehaviour());
    settings -> addOptionsBehaviour(timeUnitBehaviour);
//   settings -> addOptionsBehaviour(selectionOptionBehaviour);

    gmx::IOptionsContainer &commonOptions = options -> addGroup();
    gmx::IOptionsContainer &moduleOptions = options -> addGroup();

    settings_ -> setOptionsModuleSettings(settings);
    module_ -> initOptions(&moduleOptions, &settings_);
    settins_ -> setOptionsModuleSettings(nullptr);
// TODO: common init options
    selectionOptionsBehaviour -> initOptions(&commonOptions);
    */
}


/*
 *
 */
void
ChapRunnerModule::optionsFinished()
{
    //TODO: common_ options finished
    module_ -> optionsFinished(&settings_);
}


/*
 *
 */
int
ChapRunnerModule::run()
{
    // initialise analysis:
    //-------------------------------------------------------------------------
/*
    // common_.initTopology();
    // TODO: topology information:
    const gmx::TopologyInformation &topology;
    module_ -> initAnalysis(settings_, topology);


    // load first frame:
    //-------------------------------------------------------------------------

//    common_.initFirstFrame();
//    common_.initFrameIndexGroup();
//    module_.initAfterFirstFrame(settings_, common_.frame());

    t_pbc pbc;
    t_pbc *ppbc = settings_.hasPBC ? &pbc : NULL;
*/

}

