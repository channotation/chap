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

