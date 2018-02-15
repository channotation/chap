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

