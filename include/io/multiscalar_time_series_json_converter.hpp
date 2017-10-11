#ifndef MUTLISCALAR_TIME_SERIES_JSON_CONVERTER_HPP
#define MUTLISCALAR_TIME_SERIES_JSON_CONVERTER_HPP

#include "rapidjson/allocators.h"
#include "rapidjson/document.h"

#include "aggregation/multiscalar_time_series.hpp"


/*!
 * \brief Converts MultiscalarTimeSeries to JSON object.
 */
class MultiscalarTimeSeriesJsonConverter
{
    public:

        // conversion functionality:
        rapidjson::Value convert(
                const MultiscalarTimeSeries &timeSeries,
                rapidjson::Document::AllocatorType &alloc);

    private:

};

#endif

