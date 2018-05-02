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


#include <cmath>

#include "gromacs/analysisdata/dataframe.h"

#include "external/rapidjson/stringbuffer.h"
#include "external/rapidjson/writer.h"

#include "io/analysis_data_json_frame_exporter.hpp"


/*!
 * Returns flag indicating what types of data this module can handle.
 */
int
AnalysisDataJsonFrameExporter::flags() const
{
    return efAllowMultipoint |
           efAllowMulticolumn |
           efAllowMissing |
           efAllowMultipleDataSets;
}


/*!
 * Currently, this will open a file to which JSON data will be written. If this
 * file already exists, its content will be deleted, otherwise the file will be
 * created empty. The file stream is closed before the end of this function and
 * will be reopened for each individual frame.
 */
void
AnalysisDataJsonFrameExporter::dataStarted(
        gmx::AbstractAnalysisData* /* data */)
{
    // open file and overwrite if it already exists:
    file_.open(fileName_.c_str(), std::fstream::out);

    // TODO: header?

    // close file:
    file_.close();
}


/*!
 * This function resets the internal JSON document, adds a time stamp and frame
 * number to it, and prepares the objects and arrays that will hold all data.
 * It carries out no file system operations, which are handled by 
 * frameFinished() only.
 */
void
AnalysisDataJsonFrameExporter::frameStarted(
        const gmx::AnalysisDataFrameHeader &frame)
{   
    // calling setObject will call destructor and deallocate data:
    json_.SetObject();
    rapidjson::Document::AllocatorType& allocator = json_.GetAllocator();

    // add frame number and time stamp:
    int i = frame.index();
    real t = frame.x();
    json_.AddMember("i", i, allocator);
    json_.AddMember("t", t, allocator);

    // add object for each data set:
    for(auto it = dataSetNames_.begin(); it != dataSetNames_.end(); it++)
    {
        // create an empty dataset object to be filled when points are added:
        rapidjson::Value dataSet;
        dataSet.SetObject();

        // loop over column names and add arrays for each column:
        int colIdx = std::distance(dataSetNames_.begin(), it);
        for(auto colName : columnNames_[colIdx])
        {
            // prepare array for column:
            rapidjson::Value column;
            column.SetArray();

            // add to dataset object:
            rapidjson::Value columnName(colName, allocator);
            dataSet.AddMember(columnName, column, allocator);
        }

        // add dataset to document:
        rapidjson::Value dataSetName(*it, allocator);
        json_.AddMember(dataSetName, dataSet, allocator);
    }
}


/*!
 * Will look for the correct data frame and for this point set in the JSON
 * document prepared by frameStarted() and adds data to all column arrays. This
 * function does not perform any file system operations, which are handled by
 * frameFinished() only.
 *
 * \todo Handle errors where data set or column names are missing.
 */
void
AnalysisDataJsonFrameExporter::pointsAdded(
        const gmx::AnalysisDataPointSetRef &points)
{
    // create an allocator:
    rapidjson::Document::AllocatorType& allocator = json_.GetAllocator();

    // obtain name of data set:
    std::string dataSetName = dataSetNames_.at(points.dataSetIndex());    

    // loop over all columns:
    for(size_t i = 0; i < points.values().size(); i++)
    {
        // obtain name of column:
        std::string columnName = columnNames_.at(points.dataSetIndex()).at(i);

        // sanity check:
        if( std::isnan( points.values().at(i).value() ) )
        {
            throw std::runtime_error("Data value " + dataSetName + "/" + columnName + " is NaN and can not be "
                                     "written to JSON file.");
        }

        // add value to column array:
        rapidjson::Value val( points.values().at(i).value() );
        json_[dataSetName][columnName].PushBack(val, allocator);
    }   
}


/*!
 * Opens the file to which the JSON data will be written. The file is opened
 * in append mode after having been created/cleared by dataStarted(). The 
 * JSON document created in frameStarted() is stringified and the resulting 
 * string is streamed to the file. A new line character is also added. 
 *
 * This function handles all file system operations (i.e. opening and closing 
 * of the file and writing a string to it) and does not manipulate the JSON
 * document prepared by startFrame() and pointsAdded(). Spliiting the 
 * functionality in this way should make clean error handling possible.
 *
 * \todo Check validity of JSON document before writing.
 * \todo Check that file write was successful.
 */
void
AnalysisDataJsonFrameExporter::frameFinished(
        const gmx::AnalysisDataFrameHeader& /*frame*/)
{
    // open output file separately for each frame:
    file_.open(fileName_.c_str(), std::fstream::app);

    // create buffer writer string writer:
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    json_.Accept(writer);

    // stringify JSON document and add it to file as new line:
    std::string jsonLine(buffer.GetString(), buffer.GetSize());
    file_<<jsonLine<<std::endl;

    // TODO: should probably check if write was successful

    // close output file:
    file_.close();
}


/*!
 * Currently this does nothing and is implemented only because this is a pure
 * virtual function of the base class.
 */
void
AnalysisDataJsonFrameExporter::dataFinished()
{

}


/*!
 * Sets the name of the file to which the data will be exported.
 */
void
AnalysisDataJsonFrameExporter::setFileName(
        const std::string &fileName)
{
    fileName_ = fileName;
}


/*!
 * Setter function for data set names. Input vector should have as many 
 * elements as the number of data sets to be handled by the exporter.
 */
void
AnalysisDataJsonFrameExporter::setDataSetNames(
        const std::vector<std::string> &dataSetNames)
{
    dataSetNames_ = dataSetNames;
}


/*!
 * Setter function for column names. Input is a vector of vectors, where the 
 * outer vector should have as many elements as the number of datasets and
 * the inner vector should have as many elements as the number of columns in 
 * the respective data set.
 */
void
AnalysisDataJsonFrameExporter::setColumnNames(
        const std::vector<std::vector<std::string>> &columnNames)
{
    columnNames_ = columnNames;
}

