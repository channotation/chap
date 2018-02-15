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


#ifndef SCALAR_TIME_SERIES_HPP
#define SCALAR_TIME_SERIES_HPP

#include <map>
#include <string>

#include <gromacs/utility/real.h>


/*!
 * \brief Data container for scalar time series.
 */
class ScalarTimeSeries
{
    friend class MultiscalarTimeSeries;

    public:

        // constructor:
        ScalarTimeSeries(const std::string variableName);

        // add data point to time series:
        void addDataPoint(const real timeStamp, const real scalarValue);

        // set name of scalar variable:
        void setVariableName(const std::string variableName);

        // get name of scalar variable:
        std::string variableName() const;

    private:
    
        // name of the scalar variable:
        std::string variableName_;

        // internal data in the (time stamp, value) tuple format:
        std::map<real, real> timeSeriesData_;
};


#endif

