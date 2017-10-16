#ifndef SPLINE_CURVE_1D_JSON_CONVERTER_HPP
#define SPLINE_CURVE_1D_JSON_CONVERTER_HPP

#include "rapidjson/document.h"

#include "geometry/spline_curve_1D.hpp"


/*
 *
 */
class SplineCurve1DJsonConverter
{
    public:
        
        static SplineCurve1D fromJson(
                rapidjson::Value &val,
                unsigned int degree);

    private:

};

#endif

