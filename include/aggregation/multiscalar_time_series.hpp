#ifndef MULTISCALAR_TIME_SERIES_HPP
#define MULTISCALAR_TIME_SERIES_HPP

#include <map>
#include <vector>
#include <string>

#include "aggregation/scalar_time_series.hpp"


/*!
 * 
 */
class MultiscalarTimeSeries
{
    friend class MultiscalarTimeSeriesJsonConverter;

    public:

        // add a time series to the container:
        void addScalarTimeSeries(const ScalarTimeSeries &timeSeries);

    private:

        // internal storage of time series:
        std::map<std::string, std::vector<real>> timeSeriesData_;
};

#endif

