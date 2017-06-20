#include <cstdio>
#include <iostream>

#include <gromacs/utility/programcontext.h>
#include <gromacs/analysisdata/abstractdata.h>
#include <gromacs/analysisdata/dataframe.h>

#include "rapidjson/filewritestream.h"
#include "rapidjson/writer.h"

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

    // sanity checks:
    if( dataSetNames_.size() != data -> dataSetCount() )
    {
        std::cerr<<"ERROR: Number of data set names must equal data set count!"<<std::endl;
        std::abort();
    }
    if( columnNames_.size() != data -> dataSetCount() )
    {
        std::cerr<<"ERROR: Need to provide column names for each data set!"<<std::endl;
        std::abort();
    }

    // create a named object for each dataset:
    for(auto it = dataSetNames_.begin(); it != dataSetNames_.end(); it++)
    {
        // name of dataset as rapidjson value:
        rapidjson::Value dataSetName(*it, allocator);

        // JSON object for dataset:
        rapidjson::Value dataSet;
        dataSet.SetArray();

        // add to json document:
        json_.AddMember(dataSetName, dataSet, allocator);

        // create an array internal to this member:
        rapidjson::Value array(rapidjson::kArrayType);
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
    // sanity check:
    if( point.columnCount() != columnNames_[point.dataSetIndex()].size() )
    {
        std::cerr<<"ERROR: Number of column names must equal column count!"<<std::endl;
        std::abort();
    }

    std::vector<std::string> colNames = columnNames_[point.dataSetIndex()];

    // get an allocator:
    rapidjson::Document::AllocatorType& allocator = json_.GetAllocator();

    //
    rapidjson::Value val;
    val.SetObject();
    for(size_t i = 0; i < point.columnCount(); i++)
    {
        point.values()[i].value();
        rapidjson::Value colName(colNames[i], allocator);
        val.AddMember(colName, point.values()[i].value(), allocator);
    }

    // add record to array:    
    json_[dataSetNames_[point.dataSetIndex()]].PushBack(val, allocator);
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
    rapidjson::Writer<rapidjson::FileWriteStream> writer(os);
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


/*
 *
 */
void
AnalysisDataJsonExporter::setColumnNames(std::vector<std::vector<std::string>> columnNames)
{
    columnNames_ = columnNames;
}

