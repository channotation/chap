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


#ifndef ANALYSIS_DATA_JSON_FRAME_EXPORTER
#define ANALYSIS_DATA_JSON_FRAME_EXPORTER

#include <fstream>
#include <memory>
#include <string>

#include "gromacs/analysisdata/datamodule.h"

#include "external/rapidjson/document.h"


/*!
 * \brief This class implements the export of analysis data to a JSON file in a 
 * per-frame fashion.
 *
 * AnalysisDataJsonFrameExporter implements an AnlysisDataModuleSerial, i.e.
 * it processes data arriving from all frames in their original order rather
 * than in the order of arrival (in a parallel setting frames may be processed
 * out of order). For each frame it will create a JSON object containing all
 * data collected for that frame as well as a time and frame number stamp. 
 * These JSON objects are written to a file such that each frame object is
 * a separate row, effectively creating a file in newline delimited JSON
 * format.
 *
 * While this format is less intuitive than a tabular format, it has several
 * advantages, namely:
 *
 * - Each frame can be handled individually and there is no need to keep
 *   data for multiple frames in memory. This should allow handling of
 *   very long trajectories even if the spline curves have many control 
 *   points e.g. a small probe step was chosen.
 * - The resulting file should be comparably small, as record labels are
 *   not repeated many times and each data set column can be treated as a
 *   JSON array.
 * - The newline delimited JSON format allows streaming, i.e. the 
 *   resulting file can be read line by line (where each line is a valid
 *   JSON document), and features of interest can be extracted without the
 *   need to keep the entire file in memory (or implement a much more 
 *   complicated parser for extracting features of interest).
 *
 * The repeated opening and closing of files may not be very efficient, but 
 * will likely not be the bottleneck of the analysis tool.
 */
class AnalysisDataJsonFrameExporter : public gmx::AnalysisDataModuleSerial
{
    public:

        // constructor and destructor:
        AnalysisDataJsonFrameExporter(){};
        ~AnalysisDataJsonFrameExporter(){};

        // interface for interacting with trajectory analysis module:
        virtual int flags() const;
        virtual void dataStarted(
                gmx::AbstractAnalysisData *data);
        virtual void frameStarted(
                const gmx::AnalysisDataFrameHeader &frame);
        virtual void pointsAdded(
                const gmx::AnalysisDataPointSetRef &points);
        virtual void frameFinished(
                const gmx::AnalysisDataFrameHeader &frame);
        virtual void dataFinished();

        // setter functions for names:
        void setFileName(
                const std::string &fileName);
        void setDataSetNames(
                const std::vector<std::string> &dataSetNames);
        void setColumnNames(
                const std::vector<std::vector<std::string>> &columnNames);
    

    private:

        // names of data sets and columns:
        std::vector<std::string> dataSetNames_;
        std::vector<std::vector<std::string>> columnNames_;

        // internal variables:
        rapidjson::Document json_;
        std::string fileName_ = "stream.json";
        std::fstream file_;
};


/*!
 * Shorthand notation for smart pointer to AnalysisDataJsonFrameExporter.
 */
typedef std::shared_ptr<AnalysisDataJsonFrameExporter> AnalysisDataJsonFrameExporterPointer;

#endif

