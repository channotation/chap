#ifndef SUMMARY_STATISTICS_JSON_CONVERTER_HPP
#define SUMMARY_STATISTICS_JSON_CONVERTER_HPP

#include "rapidjson/allocators.h"
#include "rapidjson/document.h"

#include "statistics/summary_statistics.hpp"


/*!
 * \brief Converts SummaryStatistics to JSON object.
 */
class SummaryStatisticsJsonConverter
{
    public:

        // conversion functionality:
        static rapidjson::Value convert(
                const SummaryStatistics &sumStats,
                rapidjson::Document::AllocatorType &alloc);

    private:



};

#endif

