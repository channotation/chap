#include "io/spline_curve_1D_json_converter.hpp"


/*
 *
 */
SplineCurve1D
SplineCurve1DJsonConverter::fromJson(
        rapidjson::Document &doc,
        unsigned int degree) const
{

    // sanity checks:
    if( !doc.HasMember("knots") || 
        !doc["knots"].IsArray() )
    {
        throw std::logic_error("Can not construct 1D spline curve from JSON! " 
                               "No attribute 'knots' of type array found.");
    }
    if( !doc.HasMember("ctrl") || 
        !doc["ctrl"].IsArray() )
    {
        throw std::logic_error("Can not construct 1D spline curve from JSON! " 
                               "No attribute 'ctrl' of type array found.");
    }
    if( doc["knots"].Size() != doc["ctrl"].Size() )
    {
        throw std::logic_error("Can not construct 1D spline curve from Json! "
                               "Unqueal number of knots and ctrl points.");
    }

    // get knots and control points from JSON document:
    std::vector<real> knots;
    std::vector<real> ctrlPoints;
    for(size_t i = 0; i < doc["knots"].Size(); i++)
    {
        knots.push_back(doc["knots"][i].GetDouble());
        ctrlPoints.push_back(doc["ctrl"][i].GetDouble());
    }

    // add duplicate andpoint knots:
    for(unsigned int i = 0; i< degree; i++)
    {
        knots.push_back(knots.back());
        knots.insert(knots.begin(), knots.front());
    }

    // return spline curve object;
    return SplineCurve1D(degree, knots, ctrlPoints);
}

