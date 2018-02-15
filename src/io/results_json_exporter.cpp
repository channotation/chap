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


#include <fstream>
#include <exception>

#include "external/rapidjson/stringbuffer.h"
#include "external/rapidjson/writer.h"

#include "config/config.hpp"
#include "config/version.hpp"

#include "io/results_json_exporter.hpp"
#include "io/summary_statistics_json_converter.hpp"
#include "io/summary_statistics_vector_json_converter.hpp"


/*!
 * Constructor creates the basic structure of the output JSON and also 
 * populates the reproducibility information object. Other objects are created
 * and can be populated incrementally.
 */
ResultsJsonExporter::ResultsJsonExporter()
    : doc_()
{
    // overall document will be an object:
    doc_.SetObject();

    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // add reproducibility information:
    rapidjson::Value reproInfo = reproducibilityInformation();
    doc_.AddMember("reproducibilityInformation", reproInfo, alloc);
    
    // create a pathway summary object:
    rapidjson::Value pathwaySummary;
    pathwaySummary.SetObject();
    doc_.AddMember("pathwaySummary", pathwaySummary, alloc);

    // create a pathway profile object:
    rapidjson::Value pathwayProfile;
    pathwayProfile.SetObject();
    doc_.AddMember("pathwayProfile", pathwayProfile, alloc);

    // create a pathway time series object:
    rapidjson::Value pathwayScalarTs;
    pathwayScalarTs.SetObject();
    doc_.AddMember("pathwayScalarTimeSeries", pathwayScalarTs, alloc);
    
    // create a pathway profile time series object:
    rapidjson::Value pathwayProfileTs;
    pathwayProfileTs.SetObject();
    doc_.AddMember("pathwayProfileTimeSeries", pathwayProfileTs, alloc);
    
    // create a residue summary object:
    rapidjson::Value residueSummary;
    residueSummary.SetObject();
    doc_.AddMember("residueSummary", residueSummary, alloc);
}


/*!
 * Adds summary statistics of a named variable to the output document.
 */
void
ResultsJsonExporter::addPathwaySummary(
        std::string name,
        const SummaryStatistics &summary)
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // convert summary statistics:
    auto sumObj = SummaryStatisticsJsonConverter::convert(summary, alloc);

    // add to output document:
    doc_["pathwaySummary"].AddMember(toVal(name), sumObj, alloc);
}


/*!
 * Adds a set of support points to the output document. May only be called 
 * once.
 */
void
ResultsJsonExporter::addSupportPoints(const std::vector<real> &supportPoints)
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create JSON arrays for the various summary statistics:
    rapidjson::Value supPts(rapidjson::kArrayType);

    // loop over the profile and fill JSON arrays:
    for(auto s : supportPoints)
    {
        supPts.PushBack(s, alloc);
    }

    // add to table (as individual columns:
    doc_["pathwayProfile"].AddMember(toVal("s"), supPts, alloc);    
}


/*!
 * Adds a new profile to the output document. The various summary statistics
 * (i.e. min, max, mean, and standard deviation) are added as individual 
 * columns.
 *
 * Note that this requires that addSupportPoints() has already been called and
 * that the number of data points in the profile is equal to the number of 
 * support points.
 */
void
ResultsJsonExporter::addPathwayProfile(
        std::string name,
        const std::vector<SummaryStatistics> &profile)
{
    // sanity checks:
    if( !doc_["pathwayProfile"].HasMember("s") )
    {
        throw std::logic_error("Can not add profile to JSON document before "
                               "support points have been added.");
    }
    if( profile.size() != doc_["pathwayProfile"]["s"].Size() )
    {
        throw std::logic_error("Number of data points in profile must equal "
                               "number of suppoert points.");
    }

    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create JSON arrays for the various summary statistics:
    rapidjson::Value min(rapidjson::kArrayType);
    rapidjson::Value max(rapidjson::kArrayType);
    rapidjson::Value mean(rapidjson::kArrayType);
    rapidjson::Value sd(rapidjson::kArrayType);

    // loop over the profile and fill JSON arrays:
    for(auto p : profile)
    {
        min.PushBack(p.min(), alloc);
        max.PushBack(p.max(), alloc);
        mean.PushBack(p.mean(), alloc);
        sd.PushBack(p.sd(), alloc);
    }

    // add to table (as individual columns:
    doc_["pathwayProfile"].AddMember(toVal(name + "Min"), min, alloc);
    doc_["pathwayProfile"].AddMember(toVal(name + "Max"), max, alloc);
    doc_["pathwayProfile"].AddMember(toVal(name + "Mean"), mean, alloc);
    doc_["pathwayProfile"].AddMember(toVal(name + "Sd"), sd, alloc);
}


/*!
 * Adds common time stamps for all scalar time series to output document.
 */
void
ResultsJsonExporter::addTimeStamps(
        const std::vector<real> &timeStamps)
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create JSON arrays for the various summary statistics:
    rapidjson::Value tmStps(rapidjson::kArrayType);

    // loop over the profile and fill JSON arrays:
    for(auto t : timeStamps)
    {
        tmStps.PushBack(t, alloc);
    }

    // add to table (as individual columns:
    doc_["pathwayScalarTimeSeries"].AddMember(toVal("t"), tmStps, alloc); 
}


/*!
 * Adds a named scalar time series to the output document.
 */
void
ResultsJsonExporter::addPathwayScalarTimeSeries(
        std::string name,
        const std::vector<real> &timeSeries)
{
    // sanity checks:
    if( !doc_["pathwayScalarTimeSeries"].HasMember("t") )
    {
        throw std::logic_error("Can not add time series data before adding "
                               "time stamps.");
    }
    if( timeSeries.size() != doc_["pathwayScalarTimeSeries"]["t"].Size() )
    {
        throw std::logic_error("Time series must have as many data points "
                               "as there are time stamp values.");
    }

    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create JSON arrays for the various summary statistics:
    rapidjson::Value ts(rapidjson::kArrayType);

    // loop over the profile and fill JSON arrays:
    for(auto val : timeSeries)
    {
        ts.PushBack(val, alloc);
    }

    // add to table (as individual columns):
    doc_["pathwayScalarTimeSeries"].AddMember(toVal(name), ts, alloc); 

}


/*!
 * Adds temporal and spatial grid points to the output for a long-format table
 * of profile data over time and space. Should only be called once.
 */
void
ResultsJsonExporter::addPathwayGridPoints(
        const std::vector<real> &timeStamps,
        const std::vector<real> &supportPoints)
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create a JSON array to hold long format time stamps and support points:
    rapidjson::Value time(rapidjson::kArrayType);
    rapidjson::Value space(rapidjson::kArrayType);

    for(auto t : timeStamps)
    {
        for(auto s : supportPoints)
        {
            time.PushBack(t, alloc);
            space.PushBack(s, alloc);
        }
    }
    
    // add both to output document:
    doc_["pathwayProfileTimeSeries"].AddMember(toVal("t"), time, alloc);
    doc_["pathwayProfileTimeSeries"].AddMember(toVal("s"), space, alloc);
}


/*!
 * Adds a vector-valued time series to the output. Requires that 
 * addPathwayGrid() has been called before and checks that the number of data
 * points in the time series is equal to the number of grid points.
 */
void
ResultsJsonExporter::addPathwayProfileTimeSeries(
        std::string name,
        const std::vector<std::vector<real>> &timeSeries)
{
    // sanity checks:
    if( !doc_["pathwayProfileTimeSeries"].HasMember("t") ||
        !doc_["pathwayProfileTimeSeries"].HasMember("s") )
    {
        throw std::logic_error("Can not at profile time series data before "
                               "setting space time grid.");
    }
    size_t numDataPoints = timeSeries.size() * timeSeries.back().size();
    if( numDataPoints != doc_["pathwayProfileTimeSeries"]["t"].Size() )
    {
        throw std::logic_error("Time series must have as many data points "
                               "as grid points.");
    }

    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create a JSON array to hold time series values as linear array:
    rapidjson::Value ts(rapidjson::kArrayType);

    // lop over time points:
    for(auto p : timeSeries)
    {
        // loop over spatial support points:
        for(auto val : p)
        {
            ts.PushBack(val, alloc);
        }
    }

    // add to long-format table:
    doc_["pathwayProfileTimeSeries"].AddMember(toVal(name), ts, alloc);
}


/*!
 * Adds time-constant residue information (residue ID, name, chain, and 
 * hydrophobicity) to output document. Can only be called once.
 */
void
ResultsJsonExporter::addResidueInformation(
        const std::vector<int> &resId,
        const ResidueInformationProvider &resInf)
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // prepare JSON arrays for time-constant residue information:
    rapidjson::Value id(rapidjson::kArrayType);
    rapidjson::Value name(rapidjson::kArrayType);
    rapidjson::Value chain(rapidjson::kArrayType);
    rapidjson::Value hydrophobicity(rapidjson::kArrayType);

    // loop over residues:
    for(auto i : resId)
    {
        // add residue information to JSON arrays:
        id.PushBack(i, alloc);
        name.PushBack(toVal(resInf.name(i)), alloc);
        chain.PushBack(toVal(resInf.chain(i)), alloc);
        hydrophobicity.PushBack(resInf.hydrophobicity(i), alloc);
    }

    // add residue information to output document:
    doc_["residueSummary"].AddMember("id", id, alloc);
    doc_["residueSummary"].AddMember("name", name, alloc);
    doc_["residueSummary"].AddMember("chain", chain, alloc);
    doc_["residueSummary"].AddMember("hydrophobicity", hydrophobicity, alloc);
}


/*!
 * Adds summary statistics of a time-dependent residue property to the output
 * document. Requires that addResidueInformation() has been called beforehand 
 * and that the number of data points equals the number of residues.
 */
void
ResultsJsonExporter::addResidueSummary(
        std::string name,
        const std::vector<SummaryStatistics> &resSummary)
{
    // sanity checks:
    if( !doc_["residueSummary"].HasMember("id") )
    {
        throw std::logic_error("Can not add summary statistics to residue "
                               "summary before residue information has been "
                               "added.");
    }
    if( resSummary.size() != doc_["residueSummary"]["id"].Size() )
    {
        throw std::logic_error("Number of data points in summary statistics "
                               "vector must equal number residues.");
    }

    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // convert the summary statistics to JSON format:
    auto v = SummaryStatisticsVectorJsonConverter::convert(resSummary, alloc);

    // add to output document:
    doc_["residueSummary"].AddMember(toVal(name), v, alloc);
}


/*!
 * Writes the JSON document to a file of the given name.
 */
void
ResultsJsonExporter::write(std::string filename)
{
    // stringify output document:
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    doc_.Accept(writer);
    std::string outLine(buffer.GetString(), buffer.GetSize());
    
    // open outgoing file stream:
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);

    // write JSON output to file:
    file<<outLine<<std::endl;

    // close out file stream:
    file.close();
}


/*!
 * Helper function that converts a standard string into a rapidjson value to 
 * be used as e.g. member name.
 */
rapidjson::Value
ResultsJsonExporter::toVal(const std::string &str)
{
    return rapidjson::Value(str.c_str(), doc_.GetAllocator());
}


/*!
 * Returns a JSON object containing the CHAP version number and call string.
 */
rapidjson::Value
ResultsJsonExporter::reproducibilityInformation()
{
    // obtain an allocator:
    rapidjson::Document::AllocatorType &alloc = doc_.GetAllocator();

    // create an object to contain reproducibility information:
    rapidjson::Value reproInfo;
    reproInfo.SetObject();

    // create an object to contain the version information:
    rapidjson::Value version;
    version.SetObject();
    version.AddMember(
            "string",
            chapVersionString(),
            alloc);
    version.AddMember(
            "major",
            chapVersionMajor(),
            alloc);
    version.AddMember(
            "minor",
            chapVersionMinor(),
            alloc);
    version.AddMember(
            "patch",
            chapVersionPatch(),
            alloc);
    version.AddMember(
            "gitHash",
            chapVersionGitHash(),
            alloc);

    // add version information and call string to reproducibility info:
    reproInfo.AddMember(
            "version",
            version,
            alloc);
    reproInfo.AddMember(
            "commandLine",
            chapCommandLine(),
            alloc);

    // return complete reproducibility information value:
    return reproInfo;
}

