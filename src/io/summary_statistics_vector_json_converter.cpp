#include "io/summary_statistics_vector_json_converter.hpp"

#include <iostream>

/*!
 * Converts given vector of summary statistics into JSON object.
 */
rapidjson::Value
SummaryStatisticsVectorJsonConverter::convert(
        const std::vector<SummaryStatistics> &sumStats,
        rapidjson::Document::AllocatorType &alloc)
{
    rapidjson::Value sumStatsObject;
    sumStatsObject.SetObject();

    // build arrays for each summary statistic:
    rapidjson::Value min(rapidjson::kArrayType);
    rapidjson::Value max(rapidjson::kArrayType);
    rapidjson::Value mean(rapidjson::kArrayType);
    rapidjson::Value sd(rapidjson::kArrayType);
    rapidjson::Value var(rapidjson::kArrayType);
    for(auto sumStat : sumStats)
    {
        min.PushBack(sumStat.min(), alloc);
        max.PushBack(sumStat.max(), alloc);
        mean.PushBack(sumStat.mean(), alloc);
        sd.PushBack(sumStat.sd(), alloc);
        var.PushBack(sumStat.var(), alloc);
    }

    // add arrays to JSON object and return:
    sumStatsObject.AddMember("min", min, alloc);
    sumStatsObject.AddMember("max", max, alloc);
    sumStatsObject.AddMember("mean", mean, alloc);
    sumStatsObject.AddMember("sd", sd, alloc);
    sumStatsObject.AddMember("var", var, alloc);
    return sumStatsObject;
}

