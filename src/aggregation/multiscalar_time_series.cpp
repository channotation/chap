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

