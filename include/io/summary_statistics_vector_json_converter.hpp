#ifndef SUMMARY_STATISTICS_VECTOR_JSON_CONVERTER_HPP
#define SUMMARY_STATISTICS_VECTOR_JSON_CONVERTER_HPP

#include <vector>

#include "rapidjson/allocators.h"
#include "rapidjson/document.h"

#include "statistics/summary_statistics.hpp"


/*!
 * \brief Converts vector of SummaryStatistics to JSON object.
 */
class SummaryStatisticsVectorJsonConverter
{
    public:

        // conversion functionality:
        static rapidjson::Value convert(
                const std::vector<SummaryStatistics> &sumStats,
                rapidjson::Document::AllocatorType &alloc);
};

#endif

