#include "aggregation/scalar_time_series.hpp"


/*!
 * Constructor automatically sets the name of the scalar variable. 
 */
ScalarTimeSeries::ScalarTimeSeries(const std::string variableName)
    : variableName_(variableName)
{

}


/*!
 * Adds a data point to the time series. 
 *
 * \throws Logic error if time stamp already exists in time series.
 */
void
ScalarTimeSeries::addDataPoint(const real timeStamp, const real scalarValue)
{
    // check if value is already present:
    if( timeSeriesData_.find(timeStamp) != timeSeriesData_.end() )
    {
        throw std::logic_error("Time stamps must be unique in time series.");
    }
    else
    {
        // add time stamp value pair
        timeSeriesData_[timeStamp] = scalarValue;
    }
}


/*!
 * Sets the name of the scalar variable in the time series.
 */
void
ScalarTimeSeries::setVariableName(const std::string variableName)
{
    variableName_ = variableName;
}


/*!
 * Returns name of scalar variable.
 */
std::string
ScalarTimeSeries::variableName() const
{
    return variableName_;
}

