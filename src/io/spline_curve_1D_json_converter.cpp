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

