#ifndef SPLINE_CURVE_1D_JSON_CONVERTER_HPP
#define SPLINE_CURVE_1D_JSON_CONVERTER_HPP

#include "external/rapidjson/document.h"

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Conversion between a SplineCurve1D object and its JSON serialisation.
 *
 * \todo Currently onlz conversion from JSON to object implemented.
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

