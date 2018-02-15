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


#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h> 

#include "io/molecular_path_obj_exporter.hpp"

#include "geometry/cubic_spline_interp_3D.hpp"


/*!
 * Constructor initialises the \f$ s \f$ and \f$ \phi \f$ coordinates of the 
 * grid.
 */
RegularVertexGrid::RegularVertexGrid(
        std::vector<real> s,
        std::vector<real> phi)
    : s_(s)
    , phi_(phi)
{
    
}


/*!
 * Adds a vertex at the given coordinates for the given property.
 */
void
RegularVertexGrid::addVertex(
        size_t i, 
        size_t j,
        std::string p,
        gmx::RVec vertex, 
        real weight)
{
    // TODO: this situation should really be handled by a NaN colour
    if( std::isnan(weight) )
    {
        weight = 0.5;
    }

    p_.insert(p);
    std::tuple<size_t, size_t, std::string> key(i, j, p);
    vertices_[key] = vertex;
    weights_[key] = weight;
}


/*!
 * Returns a vector of all vertices for a given property.
 */
std::vector<gmx::RVec>
RegularVertexGrid::vertices(
        std::string p)
{
    // build linearly indexed vector from dual indexed map:
    std::vector<gmx::RVec> vert;
    vert.reserve(vertices_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::tuple<size_t, size_t, std::string> key(i, j, p);
            if( vertices_.find(key) != vertices_.end() )
            {
                vert.push_back(vertices_[key]); 
            }
            else
            {
                throw std::logic_error("Invalid vertex reference encountered.");
            }
        }
    }

    return vert;
}


/*!
 * Returns vector of vertex normals for a given property.
 */
std::vector<gmx::RVec>
RegularVertexGrid::normals(
        std::string p)
{
    std::vector<gmx::RVec> norm;
    norm.reserve(normals_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::tuple<size_t, size_t, std::string> key(i, j, p);
            if( normals_.find(key) != normals_.end() )
            {
                norm.push_back(normals_[key]); 
            }
            else
            {
                throw std::logic_error("Invalid vertex normal reference "
                                       "encountered.");
            }
        }
    }

    return norm;
}


/*!
 * Calculates vertex normals from triangular faces.
 */
void
RegularVertexGrid::normalsFromFaces()
{
    // array sizes for index wrap:
    int mI = s_.size();
    int mJ = phi_.size();

    for(auto p : p_)
    {
        for(int i = 0; i < s_.size(); i++)
        {
            for(int j = 0; j < phi_.size(); j++)
            {
                // index pairs for neighbouring vertices:
                std::tuple<size_t, size_t, std::string> crntKey(
                        i, 
                        j,
                        p);
                std::tuple<size_t, size_t, std::string> leftKey(
                        i, 
                        (j - 1 + mJ) % mJ,
                        p);
                std::tuple<size_t, size_t, std::string> rghtKey(
                        i, 
                        (j + 1 + mJ) % mJ,
                        p);
                std::tuple<size_t, size_t, std::string> upprKey(
                        (i + 1 + mI) % mI, 
                        j,
                        p);
                std::tuple<size_t, size_t, std::string> lowrKey((
                        i - 1 + mI) % mI, 
                        j,
                        p);
                std::tuple<size_t, size_t, std::string> dglrKey(
                        (i - 1 + mI) % mI, 
                        (j + 1 + mJ) % mJ,
                        p);
                std::tuple<size_t, size_t, std::string> dgulKey(
                        (i + 1 + mI) % mI, 
                        (j - 1 + mJ) % mJ,
                        p);

                // handle endpoints in direction along spline:
                if( i == 0 )
                {
                    lowrKey = crntKey; 
                    dglrKey = crntKey;
                }
                if( i == s_.size() - 1 )
                {
                    upprKey = crntKey;
                    dgulKey = crntKey;
                }

                // neighbouring vertices:
                gmx::RVec crntVert = vertices_.at(crntKey);
                gmx::RVec leftVert = vertices_.at(leftKey);
                gmx::RVec rghtVert = vertices_.at(rghtKey);
                gmx::RVec upprVert = vertices_.at(upprKey);
                gmx::RVec lowrVert = vertices_.at(lowrKey);
                gmx::RVec dglrVert = vertices_.at(dglrKey);
                gmx::RVec dgulVert = vertices_.at(dgulKey);

                // initialise normal as null vector:
                gmx::RVec norm(0.0, 0.0, 0.0);
                gmx::RVec sideA;
                gmx::RVec sideB;

                // North-East triangle:
                rvec_sub(rghtVert, crntVert, sideA);
                rvec_sub(upprVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);
                
                // North-North-West triangle:
                rvec_sub(upprVert, crntVert, sideA);
                rvec_sub(dgulVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);

                // West-North-West triangle:
                rvec_sub(dgulVert, crntVert, sideA);
                rvec_sub(leftVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);

                // South-West triangle:
                rvec_sub(leftVert, crntVert, sideA);
                rvec_sub(lowrVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);

                // South-South-East triangle:
                rvec_sub(lowrVert, crntVert, sideA);
                rvec_sub(dglrVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);
                
                // East-South-East triangle:
                rvec_sub(dglrVert, crntVert, sideA);
                rvec_sub(rghtVert, crntVert, sideB);
                addTriangleNorm(sideA, sideB, norm);

                // normalise normal:
                unitv(norm, norm);
                
                // add to container of normals:
                normals_[crntKey] = norm;
            }
        }
    }
}


/*!
 * Returns a vector of vertices plus a scalar weight.
 */
std::vector<std::pair<gmx::RVec, real>>
RegularVertexGrid::weightedVertices(
        std::string p)
{
    // build linearly indexed vector from dual indexed map:
    std::vector<std::pair<gmx::RVec, real>> vert;
    vert.reserve(vertices_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::tuple<size_t, size_t, std::string> key(i, j, p);
            if( vertices_.find(key) != vertices_.end() && 
                weights_.find(key) != weights_.end() )
            {
                std::pair<gmx::RVec, real> vertex(vertices_[key], weights_[key]);
                vert.push_back(vertex); 
            }
            else
            {
                throw std::logic_error("Invalid vertex reference encountered.");
            }
        }
    }

    return vert;
}


/*!
 * Returns colour scale for the given property.
 */
ColourScale
RegularVertexGrid::colourScale(std::string p)
{
    return colourScales_.at(p);
}


/*!
 * Auxiliary method for adding the contribution of two triangle sides to a
 * vertex normal.
 */
void
RegularVertexGrid::addTriangleNorm(
        const gmx::RVec &sideA,
        const gmx::RVec &sideB,
        gmx::RVec &norm)
{
    gmx::RVec tmp;
    cprod(sideA, sideB, tmp);
    rvec_add(norm, tmp, norm);
}


/*!
 * Calculates and returns vector of triangular faces for the given property.
 */
std::vector<WavefrontObjFace>
RegularVertexGrid::faces(
        std::string p)
{
    // sanity checks:
    if( phi_.size() * s_.size() * p_.size() != vertices_.size() )
    {
        throw std::logic_error("RegularVertexGrid cannot generate faces "
                               "on incomplete grid.");
    }

    if( !normals_.empty() && normals_.size() != vertices_.size() )
    {
        throw std::logic_error("Number of vertex normals does not equal "
                               "number of vertices in RegularVertexGrid.");
    }
    
    // find scalar property data range:
    real minRange = std::numeric_limits<real>::max();
    real maxRange = std::numeric_limits<real>::min();
    for(auto it = weights_.begin(); it != weights_.end(); it++)
    {
        if( it -> second < minRange )
        {
            minRange = it -> second;
        }
        if( it -> second > maxRange )
        {
            maxRange = it -> second;
        }
    }

    // prepare colour scale:
    ColourScale colScale(p);
    colScale.setRange(minRange, maxRange);
    colScale.setResolution(100);   // NOTE: limited by number of MTL materials
    colourScales_.insert(std::pair<std::string, ColourScale>(p, colScale));

    // number of vertices per property grid:
    size_t propIdx = std::distance(p_.begin(), find(p_.begin(), p_.end(), p));
    size_t vertOffset = s_.size() * phi_.size() * propIdx;

    // preallocate face vector:
    std::vector<WavefrontObjFace> faces;
    faces.reserve(phi_.size()*s_.size());

    // loop over grid:
    for(size_t i = 0; i < s_.size() - 1; i++)
    {
        for(size_t j = 0; j < phi_.size() - 1; j++)
        {
            // calculate linear indices:
            int kbl = vertOffset + i * phi_.size() + j; 
            int kbr = kbl + 1;
            int ktl = kbl + phi_.size();
            int ktr = kbr + phi_.size();

            // vertex keys:
            std::tuple<size_t, size_t, std::string> keyBl(i,     j,     p);
            std::tuple<size_t, size_t, std::string> keyBr(i,     j + 1, p);
            std::tuple<size_t, size_t, std::string> keyTl(i + 1, j,     p);
            std::tuple<size_t, size_t, std::string> keyTr(i + 1, j + 1, p);

            // face weight is average of vertex weights:
            real scalarA = weights_.at(keyBl) 
                         + weights_.at(keyTr) 
                         + weights_.at(keyTl);
            real scalarB = weights_.at(keyBl) 
                         + weights_.at(keyBr) 
                         + weights_.at(keyTr);
            scalarA /= 3.0;
            scalarB /= 3.0;

            // name of material from colour scale:
            std::string mtlNameA = colScale.scalarToColourName(scalarA); 
            std::string mtlNameB = colScale.scalarToColourName(scalarB); 

            // two faces per square:
            if( normals_.empty() )
            {
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, ktr + 1, ktl + 1},
                        mtlNameA) );
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, kbr + 1, ktr + 1},
                        mtlNameB) );
            }
            else
            {
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, ktr + 1, ktl + 1},
                        {kbl + 1, ktr + 1, ktl + 1},
                        mtlNameA) );
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, kbr + 1, ktr + 1},
                        {kbl + 1, kbr + 1, ktr + 1},
                        mtlNameB) );
            }
        }
    }

    // wrap around:
    for(size_t i = 0; i < s_.size() - 1; i++)
    {
        // calculate linear indices:
        int kbl = vertOffset + i*phi_.size() + phi_.size() - 1; 
        int kbr = vertOffset + i*phi_.size();
        int ktl = kbl + phi_.size();
        int ktr = kbr + phi_.size();

        // vertex keys:
        std::tuple<size_t, size_t, std::string> keyBl(i, phi_.size() - 1, p);
        std::tuple<size_t, size_t, std::string> keyBr(i, 0, p);
        std::tuple<size_t, size_t, std::string> keyTl(i + 1, phi_.size() - 1, p);
        std::tuple<size_t, size_t, std::string> keyTr(i + 1, 0, p);

        // face weight is average of vertex weights:
        real scalarA = weights_[keyBl] + weights_[keyTr] + weights_[keyTl];
        real scalarB = weights_[keyBl] + weights_[keyBr] + weights_[keyTr];
        scalarA /= 3.0;
        scalarB /= 3.0;

        // name of material from colour scale:
        std::string mtlNameA = colScale.scalarToColourName(scalarA); 
        std::string mtlNameB = colScale.scalarToColourName(scalarB); 

        // two faces per square:
        if( normals_.empty() )
        {
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, ktr + 1, ktl + 1},
                    mtlNameA) );
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, kbr + 1, ktr + 1},
                    mtlNameB) );
        }
        else
        {

            faces.push_back( WavefrontObjFace(
                    {kbl + 1, ktr + 1, ktl + 1},
                    {kbl + 1, ktr + 1, ktl + 1},
                    mtlNameA) );
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, kbr + 1, ktr + 1},
                    {kbl + 1, kbr + 1, ktr + 1},
                    mtlNameB) );
        }
    }

    // return face vector:
    return faces;
}



/*
 * Constructor sets default values for parameters.
 */
MolecularPathObjExporter::MolecularPathObjExporter()
    : extrapDist_(0.0)
    , gridSampleDist_(-1.0)
    , correctionThreshold_(0.1)
{
    
}


/*!
 * Sets the extrapolation distance, i.e. the distance beyond the pathway 
 * endpoints for which the surface is rendered.
 */
void
MolecularPathObjExporter::setExtrapDist(real extrapDist)
{
    extrapDist_ = extrapDist;
}


/*!
 * Sets the distance step along the along the arc of the centreline for
 * sampling surface points which are subsequently interpolated to a smooth 
 * surface.
 */
void
MolecularPathObjExporter::setGridSampleDist(real gridSampleDist)
{
    gridSampleDist_ = gridSampleDist;
}


/*!
 * Sets the threshold for identifying clashes. Is checked to lie in the 
 * interval between -1 and +1 (exclusively!), but should really be positive 
 * to identify clashes.
 */
void
MolecularPathObjExporter::setCorrectionThreshold(real correctionThreshold)
{
    // sanity check:
    if( correctionThreshold <= -1.0 + std::numeric_limits<real>::epsilon() or 
        correctionThreshold_ >= 1.0 - std::numeric_limits<real>::epsilon() )
    {
        throw std::runtime_error("Corrrection threshold for "
                                 "MolecularPathObjExporter must be between "
                                 "-1 and +1 (exclusive)!");
    }

    correctionThreshold_ = correctionThreshold;
}


/*!
 * High level driver for exporting a MolecularPath object to an OBJ and MTL
 * file.
 */
void
MolecularPathObjExporter::operator()(
        std::string fileName,
        std::string objectName,
        MolecularPath &molPath,
        std::map<std::string, ColourPalette> palettes)
{
    // define evaluation range:   
    std::pair<real, real> range(molPath.sLo() - extrapDist_,
                                molPath.sHi() + extrapDist_);

    // define resolution:
    // TODO: make this a parameter?
    int numPhi = 50;
    int numLen = std::pow(2, 8) + 1;
    std::pair<size_t, size_t> resolution(numLen, numPhi);
    
    // pathway geometry:
    auto centreLine = molPath.centreLine();
    auto pathRadius = molPath.pathRadius();

    // pathway properties:
    // (radius is added here to ensure that there is always one property)
    molPath.addScalarProperty("radius", pathRadius, false);
    auto properties = molPath.scalarProperties();   


    // Build OBJ & MTL Objects of Coloured Pore Surface
    //-------------------------------------------------------------------------

    // prepare objects:
    WavefrontObjObject obj(objectName);
    WavefrontMtlObject mtl;
  
    // generate the vertex grid:
    RegularVertexGrid grid = generateGrid(
            centreLine,
            pathRadius,
            properties,
            resolution,
            range);

    // loop over properties:
    for(auto prop : properties)
    {
        // obtain vertices, normals, and faces from grid:
        grid.normalsFromFaces();
        auto vertices = grid.weightedVertices(prop.first);
        auto vertexNormals = grid.normals(prop.first);
        auto faces = grid.faces(prop.first);

        // add faces to surface:
        WavefrontObjGroup group(prop.first);
        for(auto face : faces)
        {
            group.addFace(face);
        }

        // add to the overall OBJ object:
        obj.addVertices(vertices);
        obj.addVertexNormals(vertexNormals);
        obj.addGroup(group);

        // obtain colour scale for this property:
        auto colScale = grid.colourScale(prop.first);

        // is there a colour palatte for this property?
        if( palettes.find(prop.first) != palettes.end() )
        {
            colScale.setPalette(palettes.at(prop.first));
        }
        else if( palettes.find("default") != palettes.end() )
        {
            colScale.setPalette(palettes.at("default"));
        }
        else
        {
            throw std::runtime_error("Could not find colour palette for "
                                     "property " + prop.first + " and no "
                                     "default colour palette is available.");
        }

        // retrieve RGB colours from scale:
        auto colours = colScale.getColours();

        // create material for each colour in this scale:
        for(auto col : colours)
        {
            // create material corresponding to colour:
            WavefrontMtlMaterial material(col.first);
            material.setAmbientColour(col.second);
            material.setDiffuseColour(col.second);
            material.setSpecularColour(col.second);

            // add to material file:
            mtl.addMaterial(material);
        }
    }


    // Serialise OBJ & MTL Objects
    //-------------------------------------------------------------------------
    
    // add file extensions to base name:
    std::string objFileName = fileName + ".obj";
    std::string mtlFileName = fileName + ".mtl";

    // scale object by factor of 10 to convert nm to Ang:
    obj.scale(10.0);
    obj.calculateCog();

    // add name of material library:
    obj.setMaterialLibrary(mtlFileName);

    // create OBJ exporter and write to file:
    WavefrontObjExporter objExp;
    objExp.write(objFileName, obj);

    // create an MTL exporter and write to file:
    WavefrontMtlExporter mtlExp;
    mtlExp.write(mtlFileName, mtl);
}


/*!
 * Creates a regular vertex grid from a given centre line and radius spline.
 * This function loops over all given properties and for each property calls
 * generatePropertyGrid().
 */
RegularVertexGrid
MolecularPathObjExporter::generateGrid(
        SplineCurve3D &centreLine,
        SplineCurve1D &radius,
        std::map<std::string, std::pair<SplineCurve1D, bool>> &properties,
        std::pair<size_t, size_t> resolution,
        std::pair<real, real> range)
{
    // extract resolution:
    size_t numLen = resolution.first;
    size_t numPhi = resolution.second;

    // check if number of intervals is power of two:
    int numInt = numLen - 1;
    if( (numInt & (numInt)) == 0 && numInt > 0 )
    {
        throw std::logic_error("Number of steps along pore must be power of "
                               "two.");
    }

    // generate grid coordinates:
    std::vector<real> s;
    s.reserve(numLen);
    for(int i = 0; i < numLen; i++)
    {
        s.push_back(i*(range.second - range.first)/numLen + range.first);
    }
    std::vector<real> phi;
    phi.reserve(numPhi);
    for(size_t i = 0; i < numPhi; i++)
    {
        phi.push_back(i*2.0*M_PI/numPhi);
    }

    // generate grid from coordinates:
    RegularVertexGrid grid(s, phi);

    // loop over properties and add vertices:
    for(auto &prop : properties)
    {
        generatePropertyGrid(
                centreLine,
                radius,
                prop,
                grid);
    }

    // return the overall grid:
    return grid;
}


/*!
 * This function creates the actual vertices used in the RegularVertexGrid and
 * calculates the colour property by sampling the given spline curve at the
 * appropriate location.
 */
void
MolecularPathObjExporter::generatePropertyGrid(
        SplineCurve3D &centreLine,
        SplineCurve1D &radius,
        std::pair<std::string, std::pair<SplineCurve1D, bool>> property,
        RegularVertexGrid &grid)
{   
    // extract grid coordinates:
    std::vector<real> s;
    int num = std::floor((grid.s_.back() - grid.s_.front()) / gridSampleDist_);
    real ds = (grid.s_.back() - grid.s_.front())/(num - 1);
    for(size_t i = 0; i < num; i++)
    {
        s.push_back(grid.s_.front() + i*ds);
    }
    std::vector<real> phi = grid.phi_;
    size_t numLen = s.size();

    // sample points, radii, and tangents along molecular path:
    std::vector<gmx::RVec> centres;
    centres.reserve(s.size());
    std::vector<gmx::RVec> tangents;
    tangents.reserve(s.size());
    std::vector<real> radii;
    radii.reserve(s.size());
    for(auto eval : s)
    {
        centres.push_back( centreLine.evaluate(eval, 0) );
        gmx::RVec tv = centreLine.tangentVec(eval);
        unitv(tv, tv);
        tangents.push_back(tv);
        radii.push_back( radius.evaluate(eval, 0) );
    }

    // sample normals along molecular path:
    // (this is better then using the standard normal vector of the centre line
    // curve as it prevents "twisting" the normal vector around the curve)
    auto normals = generateNormals(tangents);
    

    // calculate sample points on pathway:
    // ------------------------------------------------------------------------

    // ring of vertices for each pathway coordinate:
    std::vector<gmx::RVec> vertRing(phi.size());
    std::map<int, std::vector<gmx::RVec>> vertexRings;
    
    // first vertex ring:
    int idxLen = 0;
    for(size_t k = 0; k < phi.size(); k++)
    {
        // rotate normal vector:
        gmx::RVec rotNormal = rotateAboutAxis(
                normals[idxLen], 
                tangents[idxLen],
                phi[k]);

        // generate vertex:
        gmx::RVec vertex = centres[idxLen];
        vertex[XX] += radii[idxLen]*rotNormal[XX];
        vertex[YY] += radii[idxLen]*rotNormal[YY];
        vertex[ZZ] += radii[idxLen]*rotNormal[ZZ];

        // add to vertex ring:
        vertRing[k] = vertex;
    }
    vertexRings[idxLen] = vertRing;

    // last vertex ring:
    idxLen = s.size() - 1;
    for(size_t k = 0; k < phi.size(); k++)
    {
        // rotate normal vector:
        gmx::RVec rotNormal = rotateAboutAxis(
                normals[idxLen], 
                tangents[idxLen],
                phi[k]);

        // generate vertex:
        gmx::RVec vertex = centres[idxLen];
        vertex[XX] += radii[idxLen]*rotNormal[XX];
        vertex[YY] += radii[idxLen]*rotNormal[YY];
        vertex[ZZ] += radii[idxLen]*rotNormal[ZZ];

        // add to vertex ring:
        vertRing[k] = vertex;
    }
    vertexRings[idxLen] = vertRing;


    // build intermediate vertex rings:    
    for(int i = 1; i <= numLen; i *= 2)
    {
        for(int j = 1; j < i; j += 2)
        {
            idxLen = j*(numLen - 1)/i;
            int idxLower = (j - 1)*(numLen - 1)/i;
            int idxUpper = (j + 1)*(numLen - 1)/i;

            // create ring of vertices:
            bool hasClashes = false;
            for(size_t k = 0; k < phi.size(); k ++)
            {
                // rotate normal vector:
                gmx::RVec rotNormal = rotateAboutAxis(
                        normals[idxLen], 
                        tangents[idxLen],
                        phi[k]);

                // generate vertex:
                gmx::RVec vertex = centres[idxLen];
                vertex[XX] += radii[idxLen]*rotNormal[XX];
                vertex[YY] += radii[idxLen]*rotNormal[YY];
                vertex[ZZ] += radii[idxLen]*rotNormal[ZZ];

                // TODO: the below goes in separate function?

                // difference vectors in neighbouring discs:
                gmx::RVec a;
                rvec_sub(vertex, centres[idxLower], a);
                unitv(a, a);
                gmx::RVec b;
                rvec_sub(vertex, centres[idxUpper], b);
                unitv(b, b);

                // cosine of angles between difference vectors:
                real cosA = iprod(a, tangents[idxLower]);
                real cosB = iprod(b, tangents[idxUpper]);

                // check overlap:
                if( cosA < correctionThreshold_ or 
                    cosB > -correctionThreshold_ )
                {
                    // set crash flag to true and terminate loop:
                    hasClashes = true;
                    break;
                }

                // add to vertex ring:
                vertRing[k] = vertex;
            }
   
            // will ignore all vertex rings with clashes:
            if( !hasClashes )
            {
                vertexRings[idxLen] = vertRing;
            }
            hasClashes = false;
        }
    }


    // interpolate 
    // ------------------------------------------------------------------------

    // interpolate support points on each equal-phi line:
    std::vector<SplineCurve3D> curves;
    for(size_t k = 0; k < phi.size(); k++)
    {
        // extract support points and parameterisation for interpolation:
        std::vector<real> param;
        std::vector<gmx::RVec> points;
        for(auto vr = vertexRings.begin(); vr != vertexRings.end(); vr++)
        {
            // pathway coordinate as curve parameter:
            param.push_back( s[vr -> first] );

            // sample points to interpolate:
            points.push_back( vr -> second[k] );
        }       

        // interpolate these points:
        CubicSplineInterp3D interp;
        curves.push_back( interp(param, points, eSplineInterpBoundaryHermite) );
    }


    // build mesh
    // ------------------------------------------------------------------------

    // sample scalar property along the path and rescale to unit interval:
    std::vector<real> prop;
    prop.reserve(grid.s_.size());
    for(size_t i = 0; i < grid.s_.size(); i++)
    {
        prop.push_back( property.second.first.evaluate(grid.s_[i], 0) );
    }
    shiftAndScale(prop, property.second.second);

    // loop over target grid coordinates and add vertices:
    for(size_t i = 0; i < grid.s_.size(); i++)
    {
        for(size_t k = 0; k < grid.phi_.size(); k++)
        {
            grid.addVertex(
                    i, 
                    k, 
                    property.first, 
                    curves[k].evaluate(grid.s_[i], 0),
                    prop[i]);
        }
    }
}


/*!
 * generates a set of normal vectors by determining a normal vector at one end
 * of the surface tube and moving it along the spline curve.
 */
std::vector<gmx::RVec>
MolecularPathObjExporter::generateNormals(
        const std::vector<gmx::RVec> &tangents)
{
    // preallocate container for normal vectors:
    std::vector<gmx::RVec> normals;
    normals.reserve(tangents.size());

    // generate initial normal:
    gmx::RVec normal = orthogonalVector(tangents[0]);
    unitv(normal, normal);
    normals.push_back(normal);
 
    // loop over tangent vectors and update normals:
    for(unsigned int i = 1; i < tangents.size(); i++)
    {
        // find axis of tangent rotation:
        gmx::RVec tangentRotAxis;
        cprod(tangents[i - 1], tangents[i], tangentRotAxis);

        // find tangent angle of rotation:
        real tangentRotCos = iprod(tangents[i], tangents[i - 1]);
        real tangentRotAxisLen = norm(tangentRotAxis);
        real tangentRotAngle = std::atan2(tangentRotAxisLen, tangentRotCos);

        // update normal by rotating it like the tangent:
        normal = rotateAboutAxis(normal, tangentRotAxis, tangentRotAngle);
        unitv(normal, normal);
        normals.push_back(normal);
    }

    // return normal vectors:
    return normals;
}


/*!
 * Returns a vector that is orthogonal to the given input vector.
 */
gmx::RVec
MolecularPathObjExporter::orthogonalVector(gmx::RVec vec)
{
    // find first nonzero element in vector:
    int idxNonZero = -1;
    for(int i = 0; i < 3; i++)
    {
        if( std::abs(vec[i]) > std::numeric_limits<real>::epsilon() )
        {
            idxNonZero = i;
            break;
        }
    }

    // sanity check:
    if( idxNonZero == -1 )
    {
        throw std::runtime_error("Can not find orthogonal to null vector.");
    }

    // find index for switching:
    int idxSwitch = (idxNonZero + 1) % 3;

    // construct non-colinear vector by element switching:
    gmx::RVec otherVec = vec;
    otherVec[idxNonZero] = vec[idxSwitch];
    otherVec[idxSwitch] = -vec[idxNonZero];

    // construct orthogonal vector via cross product:
    gmx::RVec orthVec;
    cprod(vec, otherVec, orthVec);

    // return a vector orthogonal to input vector:
    return orthVec;
}


/*!
 * Returns a vector that is rotated about the given axis by a given number of
 * degrees.
 */
gmx::RVec
MolecularPathObjExporter::rotateAboutAxis(gmx::RVec vec, 
                                          gmx::RVec axis,
                                          real angle)
{
    // evaluate trigonometric functions of rotation angle:
    const real COS = std::cos(angle);
    const real SIN = std::sin(angle);

    // construct rotation matrix:
    matrix rotMat;
    rotMat[XX][XX] = COS + axis[XX]*axis[XX]*(1.0 - COS);
    rotMat[XX][YY] = axis[XX]*axis[YY]*(1.0 - COS) - axis[ZZ]*SIN;
    rotMat[XX][ZZ] = axis[XX]*axis[ZZ]*(1.0 - COS) + axis[YY]*SIN;
    rotMat[YY][XX] = axis[YY]*axis[XX]*(1.0 - COS) + axis[ZZ]*SIN;
    rotMat[YY][YY] = COS + axis[YY]*axis[YY]*(1 - COS);
    rotMat[YY][ZZ] = axis[YY]*axis[ZZ]*(1.0 - COS) - axis[XX]*SIN;
    rotMat[ZZ][XX] = axis[ZZ]*axis[XX]*(1.0 - COS) - axis[YY]*SIN;
    rotMat[ZZ][YY] = axis[ZZ]*axis[YY]*(1.0 - COS) + axis[XX]*SIN;
    rotMat[ZZ][ZZ] = COS + axis[ZZ]*axis[ZZ]*(1.0 - COS);

    // rotate vector:
    gmx::RVec rotVec;
    mvmul(rotMat, vec, rotVec);

    // return rotated vector:
    return rotVec;
}



/*!
 * Shift ands scales all values in input vector so that they lie in the unit
 * interval.
 *
 * The divergent flag controls how the data is scaled. For a divergent scale, 
 * both positive and negative values are scaled by the same factor and then
 * shifted to the unit interval, in this case the zero of the scale is
 * precisely 0.5. For a sequential colour scale, all values are first shifted 
 * to the positive real range and then scaled to the unit interval, this does
 * obviously not preserve the zero of the original array.
 */
void
MolecularPathObjExporter::shiftAndScale(
        std::vector<real> &prop,
        bool divergent)
{
    // find data range:
    real minProp = *std::min_element(prop.begin(), prop.end()); 
    real maxProp = *std::max_element(prop.begin(), prop.end());
    
    // scale for divergent colour scale?
    if( std::fabs(maxProp - minProp) < std::numeric_limits<real>::epsilon() )
    {
        // special case of constant property value, shift to middle of scale:
        real shift = -minProp + 0.5;

        // just shift values in this case:
        std::for_each(prop.begin(), prop.end(), [shift](real &p){p += shift;});
    }
    else if( divergent == false )
    {
        // for sequential colour scale, shift data to positive real range and
        // scale by length of data interval:
        real shift = -minProp;
        real scale = 1.0/(maxProp - minProp);  
       
        // shift and scale the property array:
        std::for_each(prop.begin(), prop.end(), [shift](real &p){p += shift;});
        std::for_each(prop.begin(), prop.end(), [scale](real &p){p *= scale;});
    }
    else
    {
        // for divergent colour scale, scale both positive and negative values
        // such that data lies in [-0.5, 0.5], then shift by 0.5:
        real shift = 0.5;
        real scale = 1.0/std::max(std::fabs(minProp), std::fabs(maxProp))/2.0; 

        // scale and shift:
        std::for_each(prop.begin(), prop.end(), [scale](real &p){p *= scale;});
        std::for_each(prop.begin(), prop.end(), [shift](real &p){p += shift;});
    }
}

