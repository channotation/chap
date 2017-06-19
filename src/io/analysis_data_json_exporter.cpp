#include <cstdio>
#include <iostream>

#include <gromacs/utility/programcontext.h>
#include <gromacs/analysisdata/abstractdata.h>

#include "rapidjson/filewritestream.h"
#include "rapidjson/prettywriter.h"

#include "io/analysis_data_json_exporter.hpp"


/*
 *
 */
int
AnalysisDataJsonExporter::flags() const
{
    return efAllowMultipoint |
           efAllowMulticolumn |
           efAllowMissing |
           efAllowMultipleDataSets;
}


/*
 *
 */
bool
AnalysisDataJsonExporter::parallelDataStarted(
        gmx::AbstractAnalysisData *data,
        const gmx::AnalysisDataParallelOptions &options)
{
    // create JSOn object in document:
    json_.SetObject();

    // retrieve an allocator:
    rapidjson::Document::AllocatorType& allocator = json_.GetAllocator();


    // reproducibility information:
    //-------------------------------------------------------------------------

    // get program context:
    const gmx::IProgramContext &programContext = gmx::getProgramContext();

    // add reproducibility information to JSON document:
    rapidjson::Value reproInfo;
    reproInfo.SetObject();
    reproInfo.AddMember("commandline", 
                        std::string(programContext.commandLine()), 
                        allocator);
    json_.AddMember("reproducibility information", reproInfo, allocator);


    // build object for each data set:
    //-------------------------------------------------------------------------

    // sanity check:
    if( dataSetNames_.size() != data -> dataSetCount() )
    {
        std::cerr<<"ERROR: Number of data set names must equal data set count!"<<std::endl;
        std::abort();
    }

    // create a named object for each dataset:
    for(auto it = dataSetNames_.begin(); it != dataSetNames_.end(); it++)
    {
        // name of dataset as rapidjson value:
        rapidjson::Value dataSetName(*it, allocator);

        // JSON object for dataset:
        rapidjson::Value dataSet;
        dataSet.SetObject();

        // add to json document:
        json_.AddMember(dataSetName, dataSet, allocator);
    }


    // indicate that parallel support is enabled:
    return true;
}


/*
 *
 */
void
AnalysisDataJsonExporter::frameStarted(
        const gmx::AnalysisDataFrameHeader &frame)
{

}


/*
 *
 */
void
AnalysisDataJsonExporter::pointsAdded(
        const gmx::AnalysisDataPointSetRef &point)
{

}


/*
 *
 */
void
AnalysisDataJsonExporter::frameFinished(
        const gmx::AnalysisDataFrameHeader &frame)
{

}


/*
 *
 */
void
AnalysisDataJsonExporter::frameFinishedSerial(int index)
{

}


/*
 *
 */
void
AnalysisDataJsonExporter::dataFinished()
{
    // open output file:
    FILE* file = std::fopen("output.json", "w");

    // prepare buffer for JSOn output:
    char buffer[65536];
    rapidjson::FileWriteStream os(file, buffer, sizeof(buffer));

    // write prettyfied JSON to file:
    rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(os);
    json_.Accept(writer);

    // close output file:
    std::fclose(file);
}


/*
 *
 */
void
AnalysisDataJsonExporter::setDataSetNames(std::vector<std::string> dataSetNames)
{
    dataSetNames_ = dataSetNames;
}

