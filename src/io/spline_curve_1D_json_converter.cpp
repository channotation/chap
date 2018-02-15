#include "io/spline_curve_1D_json_converter.hpp"


/*!
 * Converts a JSON value with (unique) knots and control point arrays to a 
 * SplineCurve1D.
 *
 * \todo Degree of spline curve must currently be given explicitly.
 */
SplineCurve1D
SplineCurve1DJsonConverter::fromJson(
        rapidjson::Value &val,
        unsigned int degree)
{

    // sanity checks:
    if( !val.HasMember("knots") || 
        !val["knots"].IsArray() )
    {
        throw std::logic_error("Can not construct 1D spline curve from JSON! " 
                               "No attribute 'knots' of type array found.");
    }
    if( !val.HasMember("ctrl") || 
        !val["ctrl"].IsArray() )
    {
        throw std::logic_error("Can not construct 1D spline curve from JSON! " 
                               "No attribute 'ctrl' of type array found.");
    }
    if( val["knots"].Size() != val["ctrl"].Size() )
    {
        throw std::logic_error("Can not construct 1D spline curve from Json! "
                               "Unqueal number of knots and ctrl points.");
    }

    // get knots and control points from JSON valument:
    std::vector<real> knots;
    std::vector<real> ctrlPoints;
    for(size_t i = 0; i < val["knots"].Size(); i++)
    {
        knots.push_back(val["knots"][i].GetDouble());
        ctrlPoints.push_back(val["ctrl"][i].GetDouble());
    }

    // add duplicate endpoint knots:
    for(unsigned int i = 0; i< degree; i++)
    {
        knots.push_back(knots.back());
        knots.insert(knots.begin(), knots.front());
    }

    // return spline curve object;
    return SplineCurve1D(degree, knots, ctrlPoints);
}

