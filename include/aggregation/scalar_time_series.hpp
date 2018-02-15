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

