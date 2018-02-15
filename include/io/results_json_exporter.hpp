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


#ifndef RESULTS_JSON_EXPORTER_HPP
#define RESULTS_JSON_EXPORTER_HPP

#include <string>

#include "external/rapidjson/document.h"

#include "analysis-setup/residue_information_provider.hpp"
#include "statistics/summary_statistics.hpp"


/*!
 * \brief Container class for facilitating the export of results to a JSON file.
 */
class ResultsJsonExporter
{
    public:
       
        // constructor:
        ResultsJsonExporter();

        // interface for adding to output:
        void addPathwaySummary(
                std::string name,
                const SummaryStatistics &summary);
        void addSupportPoints(
                const std::vector<real> &supportPoints);
        void addPathwayProfile(
                std::string name,
                const std::vector<SummaryStatistics> &profile);
        void addTimeStamps(
                const std::vector<real> &timeStamps);
        void addPathwayScalarTimeSeries(
                std::string name,
                const std::vector<real> &timeSeries);
        void addPathwayGridPoints(
                const std::vector<real> &timeStamps,
                const std::vector<real> &supportPoints);
        void addPathwayProfileTimeSeries(
                std::string name,
                const std::vector<std::vector<real>> &timeSeries);
        void addResidueInformation(
                const std::vector<int> &resId,
                const ResidueInformationProvider &resInf);
        void addResidueSummary(
                std::string name,
                const std::vector<SummaryStatistics> &resSummary);

        // interface for writing to file:
        void write(std::string filename);

    private:

        // helper function to convert a string to a rapidjson value:
        inline rapidjson::Value toVal(const std::string &str);

        // function for populating the reproducibility info with values:
        rapidjson::Value reproducibilityInformation();

        // overall output document:
        rapidjson::Document doc_;
};

#endif

