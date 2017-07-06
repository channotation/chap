#include "io/summary_statistics_json_converter.hpp"


/*
 *
 *
 */
rapidjson::Value
SummaryStatisticsJsonConverter::convert(
        const SummaryStatistics &sumStats,
        rapidjson::Document::AllocatorType &alloc)
{
    // create JSON objecy to hold summary statistics:
    rapidjson::Value sumStatsObject;
    sumStatsObject.SetObject();
    
    // add different summary statistics:
    sumStatsObject.AddMember(
            "min", 
            sumStats.min(), 
            alloc);
    sumStatsObject.AddMember(
            "max", 
            sumStats.max(), 
            alloc);
    sumStatsObject.AddMember(
            "mean", 
            sumStats.mean(), 
            alloc);
    sumStatsObject.AddMember(
            "sd", 
            sumStats.sd(), 
            alloc);
    sumStatsObject.AddMember(
            "var", 
            sumStats.var(), 
            alloc); 

    // return object:
    return sumStatsObject;
}

