#include "aggregation/multiscalar_time_series.hpp"


/*!
 * Adds a time series to the container.
 */
void
MultiscalarTimeSeries::addScalarTimeSeries(const ScalarTimeSeries &timeSeries)
{
    // check if we already have time series:
    if( timeSeriesData_.empty() )
    {
        // add time stamps and scalar values to internal storage:
        std::vector<real> timeStamp;
        std::vector<real> scalarValue;
        for(auto& dataPoint : timeSeries.timeSeriesData_)
        {
            timeStamp.push_back(dataPoint.first);
            scalarValue.push_back(dataPoint.second);
        }
        timeSeriesData_["t"] = timeStamp;
        timeSeriesData_[timeSeries.variableName()] = scalarValue;
    }
    else
    {
        // check that time stamp is identical to existing ones:
        for(auto& timeStamp : timeSeriesData_["t"])
        {
            if( timeSeries.timeSeriesData_.find(timeStamp) == 
                timeSeries.timeSeriesData_.end() )
            {
                throw std::logic_error("Time stamps of new time series do "
                "not match existing time stamps.");
            }
        }

        // copy scalar values to new vector:
        std::vector<real> scalarValue;    
        for(auto& dataPoint : timeSeries.timeSeriesData_)
        {
            scalarValue.push_back(dataPoint.second);
        }

        // add only scalar values to internal storage:
        if( timeSeriesData_.find(timeSeries.variableName()) != 
            timeSeriesData_.end() )
        {
            throw std::logic_error("Multiscalar time series can not accept "
            "two scalar time series with the same variable name.");
        }
        else
        {
            timeSeriesData_[timeSeries.variableName()] = scalarValue;
        }
    }
}

