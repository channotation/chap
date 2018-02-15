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


#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <map>
#include <string>
#include <tuple>
#include <unordered_set>

#include <gtest/gtest.h>

#include "path-finding/molecular_path.hpp"
#include "io/colour.hpp"
#include "io/wavefront_mtl_io.hpp"
#include "io/wavefront_obj_io.hpp"


/*!
 * \brief Representation of a regular grid on a cylinder surface.
 *
 * This class is used as a smart container for vertices on the surface of a
 * cylinder of tube, i.e. a grid in the \f$ s \f$ and \f$ \phi \f$ space. 
 * As the grid is regular, it is know which vertices neighbour one another,
 * which can be exploited to generate triangular faces. These faces can be used
 * to subsequently generate vertex normals. 
 *
 * This is all used by MolcularPathObjExporter.
 */
class RegularVertexGrid
{
    friend class MolecularPathObjExporter;

    public:
        
        // constructor:
        RegularVertexGrid(
                std::vector<real> s,
                std::vector<real> phi);

        // interface for adding vertices to the grid:
        void addVertex(
                size_t i, 
                size_t j,
                std::string p,
                gmx::RVec vertex, 
                real weight);


        void addColourScale(
                std::string p,
                std::vector<real> prop,
                std::vector<gmx::RVec> palette);



        void interpolateMissing();
        
        // method for determining vertex normals:
        void normalsFromFaces();
    
        // getter methods:
        std::vector<gmx::RVec> vertices(
                std::string p);
        std::vector<std::pair<gmx::RVec, real>> weightedVertices(
                std::string p);
        std::vector<gmx::RVec> normals(
                std::string p);
        std::vector<WavefrontObjFace> faces(
                std::string p);
        ColourScale colourScale(
                std::string p);

    private:

        int numX;
        int numY;

        const std::vector<real> phi_;
        const std::vector<real> s_;
        std::unordered_set<std::string> p_;

        std::map<std::string, ColourScale> colourScales_;

        std::map<std::tuple<size_t, size_t, std::string>, gmx::RVec> vertices_;
        std::map<std::tuple<size_t, size_t, std::string>, real> weights_;
        std::map<std::tuple<size_t, size_t, std::string>, gmx::RVec> normals_;



        void addTriangleNorm(
                const gmx::RVec &sideA, 
                const gmx::RVec &sideB, 
                gmx::RVec &norm);
};


/*!
 * \brief Writes the surface of a MolecularPathway to an OBJ and MTL file.
 *
 * The resulting OBJ file contains different groups of faces, each representing
 * a different scalar property mapped to the pathway surface. The colour 
 * associated with this property is written to an MTL file, which is referenced
 * at the beginning of the OBJ file.
 */
class MolecularPathObjExporter
{
    friend class MolecularPathObjExporterTest;
    FRIEND_TEST(
            MolecularPathObjExporterTest, 
            MolecularPathObjExporterOrthogonalVectorTest);
    FRIEND_TEST(
            MolecularPathObjExporterTest, 
            MolecularPathObjExporterAxisRotationTest);


    public:
        
        // constructor:
        MolecularPathObjExporter();

        // setter functions:
        void setExtrapDist(real extrapDist);
        void setGridSampleDist(real gridSampleDist);
        void setCorrectionThreshold(real correctionThreshold);
        void setPermitClashes(bool permitClashes);

        // interface for exporting:
        void operator()(
                std::string fileName,
                std::string objectName,
                MolecularPath &molPath,
                std::map<std::string, ColourPalette> palettes);


    private:

        // parameters:
        real extrapDist_;
        real gridSampleDist_;
        real correctionThreshold_;

        // functions for generating the pathway surface grid:
        std::vector<gmx::RVec> generateNormals(
                const std::vector<gmx::RVec> &tangents);
        RegularVertexGrid generateGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                std::map<std::string, std::pair<SplineCurve1D, bool>> &properties,
                std::pair<size_t, size_t> resolution,
                std::pair<real, real> range);
        void generatePropertyGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                std::pair<std::string, std::pair<SplineCurve1D, bool>> property,
                RegularVertexGrid &grid);

        // auxiliary geometric functions: 
        gmx::RVec orthogonalVector(gmx::RVec vec);
        gmx::RVec rotateAboutAxis(gmx::RVec vec, gmx::RVec axis, real angle);

        // manipulate scalar property:
        void shiftAndScale(std::vector<real> &prop, bool divergent);
};


#endif

