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


#include <iomanip>
#include <sstream>

#include "io/colour.hpp"
#include "io/json_doc_importer.hpp"


/*!
 * Constructs a palette where the given values are associated with the given
 * colours.
 */
ColourPalette::ColourPalette(
        std::vector<real> values,
        std::vector<gmx::RVec> colours)
    : values_(values)
    , colours_(colours)
{
    
}


/*!
 * Copy constructor for colour palettes.
 */
ColourPalette::ColourPalette(
        const ColourPalette &other)
    : values_(other.values_)
    , colours_(other.colours_)
{
    
}


/*!
 * Creates a names colour scale.
 */
ColourScale::ColourScale(std::string name)
    : name_(name)
{
    
}


/*!
 * Sets the palette, i.e. the finite set of RGB colours which will be 
 * interpolated linearly to yield a continuous colour scale.
 */
void
ColourScale::setPalette(ColourPalette palette)
{
    palette_ = palette;
}


/*!
 * Sets the range of the variable mapped through this colour scale.
 */
void
ColourScale::setRange(real min, real max)
{
    // sanity check:
    if( min > max )
    {
        throw std::logic_error("Lower range limit must be smaller than upper "
                               "range limit.");
    }

    rangeMin_ = min;
    rangeMax_ = max;
}


/*!
 * Sets the number of distinct colours in this scale.
 */
void
ColourScale::setResolution(size_t res)
{
    numColours_ = res;

    // rebuild colour names:
    std::stringstream ss;
    colourNames_.resize(numColours_);
    for(size_t i = 0; i < numColours_; i++)
    {
        ss<<"_"<<std::setw(4)<<std::setfill('0')<<i;
        colourNames_[i] = name_ + ss.str();
        ss.str(std::string());
        ss.clear();
    }
}


/*!
 * Returns a map of named colours in the scale. Requires that setRange(), 
 * setResolution(), and setPalette() have been called previously.
 */
std::map<std::string, gmx::RVec>
ColourScale::getColours()
{
    // associate each colour in the palette with a scalar value:
    auto paletteValues = palette_.values();
    auto paletteColours = palette_.colours();

    // compute linear range of scalars:
    std::map<std::string, gmx::RVec> colours;
    for(int i = 0; i < colourNames_.size(); i++)
    {
        // scalar of this colour:
        real scalar = i*(rangeMax_ - rangeMin_)/(colourNames_.size() - 1);
        scalar += rangeMin_;

        // get interval boundaries:
        auto lb = std::lower_bound(
                paletteValues.begin(),
                paletteValues.end(),
                scalar);
        if( std::next(lb) == paletteValues.end() )
        {
            lb--;
        }
        size_t idxLo = std::distance(paletteValues.begin(), lb);
        size_t idxHi = idxLo + 1;

        // scaling factor for linear interpolation:
        real alpha = (scalar - paletteValues[idxLo]);
        alpha /= (paletteValues[idxHi] - paletteValues[idxLo]); 

        // interpolate palette colours:
        gmx::RVec colLo = paletteColours.at(idxLo);
        gmx::RVec colHi = paletteColours.at(idxHi);
        svmul(1.0 - alpha, colLo, colLo);
        svmul(alpha, colHi, colHi);
        gmx::RVec colInterp;
        rvec_add(colLo, colHi, colInterp);

        // colour from linear interpolation of palette:
        colours[colourNames_[i]] = colInterp;
    }

    return colours;
}


/*!
 * Returns the name of the colour corresponding to the given scalar.
 */
std::string
ColourScale::scalarToColourName(real scalar)
{
    // sanity checks:
    if( scalar < rangeMin_ || scalar > rangeMax_ )
    {
        throw std::logic_error("Value outside scale range encountered.");
    }

    // determine index of colour name:
    int idx = (numColours_ - 1)*(scalar - rangeMin_) / (rangeMax_ - rangeMin_);

    // return appropriate colour name:
    return colourNames_.at(idx);
}


/*!
 * Returns a map of named colour palettes from the given JSON document.
 */
std::map<std::string, ColourPalette>
ColourPaletteProvider::fromJsonDoc(const rapidjson::Document &doc)
{
    // sanity checks:
    if( !doc.IsObject() )
    {
        throw std::runtime_error("No valid JSON object provided for generation"
        " of colour palettes.");
    }
    if( !doc.HasMember("palettes") || !doc["palettes"].IsArray() )
    {
        throw std::runtime_error("JSON document provided for colour palette "
        "generation does not contain palette array.");
    }

    // extract palette array:
    const rapidjson::Value &palettes = doc["palettes"];

    // iterate over provided values:
    std::map<std::string, ColourPalette> paletteMap;
    rapidjson::Value::ConstValueIterator it;
    for(it = palettes.Begin(); it != palettes.End(); it++)
    {
        // check that required entries are present and of correct type:
        if( !(it -> HasMember("property")) || 
            !(*it)["property"].IsString() )
        {
            throw std::runtime_error("No 'property' attribute of type string "
            "in colour palette.");
        }
        if( !(it -> HasMember("palette")) || 
            !(*it)["palette"].IsObject() )
        {
            throw std::runtime_error("No 'palette' attribute of type object "
            "in colour palette.");
        }

        // extract palette:
        const rapidjson::Value &pal = (*it)["palette"];

        // check that all value and colour definitions are present:
        if( !pal.HasMember("val") || !pal["val"].IsArray() )
        {
            throw std::runtime_error("No value definition found in palette.");
        }
        if( !pal.HasMember("r") || !pal["r"].IsArray() )
        {
            throw std::runtime_error("No red colour definition found in palette.");
        }
        if( !pal.HasMember("g") || !pal["g"].IsArray() )
        {
            throw std::runtime_error("No green colour definition found in palette.");
        }
        if( !pal.HasMember("b") || !pal["b"].IsArray() )
        {
            throw std::runtime_error("No blue colour definition found in palette.");
        }

        // check array sizes:
        if( pal["val"].Size() != pal["r"].Size() )
        {
            throw std::runtime_error("Need one red value for each value.");
        }
        if( pal["val"].Size() != pal["g"].Size() )
        {
            throw std::runtime_error("Need one green value for each value.");
        }
        if( pal["val"].Size() != pal["b"].Size() )
        {
            throw std::runtime_error("Need one blue value for each value.");
        }

        // extract palette definition:
        std::vector<real> values;
        values.reserve(pal["val"].Size());
        std::vector<gmx::RVec> colours;
        colours.reserve(values.size());
        for(size_t i = 0; i < pal["val"].Size(); i++)
        {
            // extract scale value:
            values.push_back(pal["val"][i].GetDouble());

            // extract colour at value:
            gmx::RVec col(pal["r"][i].GetDouble(),
                          pal["g"][i].GetDouble(),
                          pal["b"][i].GetDouble());
            colours.push_back(col);
        }

        // add to palettes map:
        paletteMap[(*it)["property"].GetString()] = ColourPalette(values, colours);
    }

    return paletteMap;
}


/*!
 * Returns a map of named colour palettes imported from the given JSON file.
 */
std::map<std::string, ColourPalette>
ColourPaletteProvider::fromJsonFile(std::string filename)
{
    // import colour palette definitions from JSON file:
    JsonDocImporter import;
    auto doc = import(filename);

    // document parsing is handles by separate function:
    return fromJsonDoc(doc);
}

