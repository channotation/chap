#include <cstdio>
#include <iostream>

#include <gromacs/utility/programcontext.h>
#include <gromacs/analysisdata/abstractdata.h>
#include <gromacs/analysisdata/dataframe.h>

#include "rapidjson/filewritestream.h"
#include "rapidjson/writer.h"

#include "config/version.hpp"
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
    
    // get chap version:
    std::string version = chapVersionString();

    // add reproducibility information to JSON document:
    rapidjson::Value reproInfo;
    reproInfo.SetObject();
    reproInfo.AddMember("version",
                        version,
                        allocator);
    reproInfo.AddMember("commandline", 
                        std::string(programContext.commandLine()), 
                        allocator);
    json_.AddMember("reproducibility information", reproInfo, allocator);


    // parameter information:
    //-------------------------------------------------------------------------

    // create object for parameters:
    rapidjson::Value paramInfo;
    paramInfo.SetObject();

    // loop over all known parameters:
    for(auto it = parameterMap_.begin(); it != parameterMap_.end(); it++)
    {       
        rapidjson::Value key(it -> first, allocator);
        paramInfo.AddMember(key, it -> second, allocator);
    }

    // add parameter object to document:
    json_.AddMember("parameters", paramInfo, allocator);

    // residue information:
    //-------------------------------------------------------------------------

    // add all residue ID/name pairs to array:
    rapidjson::Value resNames(rapidjson::kArrayType);
    for(auto it = residueNames_.begin(); it != residueNames_.end(); it++)
    {
        rapidjson::Value res;
        res.SetObject();
        res.AddMember("res.id", it -> first, allocator);
        res.AddMember("res.name", it -> second, allocator);
        resNames.PushBack(res, allocator);
    }
    
    // add array to JSON document:
    json_.AddMember("residue.names", resNames, allocator);


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

    // create record for this point and add time stamp:
    rapidjson::Value val;
    val.SetObject();
    val.AddMember("t", point.x(), allocator);

    // loop over columns in this dataset and add values:
    for(size_t i = 0; i < point.columnCount(); i++)
    {
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
    // sanity check:
    if( fileName_.empty() )
    {
        std::cerr<<"ERROR: Output file name not given."<<std::endl;
        std::abort();
    }

    // open output file:
    FILE* file = std::fopen(fileName_.c_str(), "w");

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


/*
 *
 */
void
AnalysisDataJsonExporter::setResidueNames(std::unordered_map<int, std::string> resNames)
{
    residueNames_ = resNames;
}


/*
 * Setting function for output file name.
 */
void
AnalysisDataJsonExporter::setFileName(std::string fileName)
{
    fileName_ = fileName;
}


/*!
 * Adds a parameter key value pair that will be written to the JSON output. 
 * Needs to be called before dataStarted() or parallelDataStarted(), otherwise
 * the parameter will be ignored.
 *
 * TODO: make it so that this adds to the JSON doc directly!
 */
void
AnalysisDataJsonExporter::addParameter(std::string name, real value)
{
    parameterMap_[name] = std::to_string(value);
}


