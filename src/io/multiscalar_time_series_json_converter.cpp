#include "io/multiscalar_time_series_json_converter.hpp"


/*!
 * Converts a given MultiscalarTimeSeries object to a JSON object.
 */
rapidjson::Value
MultiscalarTimeSeriesJsonConverter::convert(
        const MultiscalarTimeSeries &timeSeries,
        rapidjson::Document::AllocatorType &alloc)
{
    // prepare JSON object for time series object:
    rapidjson::Value timeSeriesJson;
    timeSeriesJson.SetObject();

    // iterate over columns in data frame:
    for(auto& column : timeSeries.timeSeriesData_)
    {
        // create JSON array for data:
        rapidjson::Value jsonColumn(rapidjson::kArrayType);

        // copy column elements to JSON array:
        for(real value : column.second)
        {
            jsonColumn.PushBack(value, alloc);
        }

        std::string columnName("test");

        // add column to JSON object:
        timeSeriesJson.AddMember(
                rapidjson::StringRef(column.first.c_str()), 
                jsonColumn, 
                alloc);
    }
        
    // return JSON object:
    return timeSeriesJson;
}

