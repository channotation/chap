#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <map>
#include <string>
#include <tuple>
#include <unordered_set>

#include <gtest/gtest.h>

#include "path-finding/molecular_path.hpp"
#include "io/wavefront_mtl_io.hpp"
#include "io/wavefront_obj_io.hpp"


/*
 *
 */
class RegularVertexGrid
{
    friend class MolecularPathObjExporter;

    public:

        RegularVertexGrid(
                std::vector<real> s,
                std::vector<real> phi);

        void addVertex(
                size_t i, 
                size_t j,
                std::string p,
                gmx::RVec vertex, 
                real weight);
        void addVertexNormal(
                size_t i, 
                size_t j,
                std::string p,
                gmx::RVec normal);

        void interpolateMissing();
        void normalsFromFaces();
    
        std::vector<gmx::RVec> vertices(
                std::string p);
        std::vector<std::pair<gmx::RVec, real>> weightedVertices(
                std::string p);
        std::vector<gmx::RVec> normals(
                std::string p);
        std::vector<WavefrontObjFace> faces(
                std::string p);

    private:

        int numX;
        int numY;

        const std::vector<real> phi_;
        const std::vector<real> s_;
        std::unordered_set<std::string> p_;


        std::map<std::tuple<size_t, size_t, std::string>, gmx::RVec> vertices_;
        std::map<std::tuple<size_t, size_t, std::string>, real> weights_;
        std::map<std::tuple<size_t, size_t, std::string>, gmx::RVec> normals_;



        void addTriangleNorm(
                const gmx::RVec &sideA, 
                const gmx::RVec &sideB, 
                gmx::RVec &norm);
};


/*
 *
 */
class MolecularPathObjExporter
{
    friend class MolecularPathObjExporterTest;
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterOrthogonalVectorTest);
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterAxisRotationTest);


    public:
        
        // constructor:
        MolecularPathObjExporter();

        // interface for exporting:
        void operator()(std::string fileName,
                        std::string objectName,
                        MolecularPath &molPath);


    private:

        // 
        inline int numPlanarVertices(real &d, real &r);

        // 
        std::vector<gmx::RVec> generateNormals(
                const std::vector<gmx::RVec> &tangents);
        RegularVertexGrid generateGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                std::map<std::string, SplineCurve1D> &properties,
                std::pair<size_t, size_t> resolution,
                std::pair<real, real> range);
        void generatePropertyGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                std::pair<std::string, SplineCurve1D> property,
                RegularVertexGrid &grid);

        // 
        inline std::pair<std::vector<gmx::RVec>, std::vector<gmx::RVec>> vertexRing(
                gmx::RVec base,
                gmx::RVec tangent,
                gmx::RVec normal,
                real radius,
                real angleIncrement,
                size_t nIncrements);

        // 
        gmx::RVec orthogonalVector(gmx::RVec vec);
        gmx::RVec rotateAboutAxis(gmx::RVec vec, gmx::RVec axis, real angle);

        real cosAngle(const gmx::RVec &vecA, const gmx::RVec &vecB);

        void shiftAndScale(std::vector<real> &prop);
};


#endif

