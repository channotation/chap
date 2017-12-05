#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <map>
#include <string>
#include <tuple>
#include <unordered_set>

#include <gtest/gtest.h>

#include "io/wavefront_obj_io.hpp"
#include "path-finding/molecular_path.hpp"


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
                size_t p,
                gmx::RVec vertex, 
                real weight);
        void addVertexNormal(
                size_t i, 
                size_t j,
                size_t p,
                gmx::RVec normal);

        void interpolateMissing();
        void normalsFromFaces();
    
        std::vector<gmx::RVec> vertices(
                size_t p);
        std::vector<std::pair<gmx::RVec, real>> weightedVertices(
                size_t p);
        std::vector<gmx::RVec> normals(
                size_t p);
        std::vector<WavefrontObjFace> faces(
                size_t p);

    private:

        int numX;
        int numY;

        const std::vector<real> phi_;
        const std::vector<real> s_;
        std::unordered_set<size_t> p_;


        std::map<std::tuple<size_t, size_t, size_t>, gmx::RVec> vertices_;
        std::map<std::tuple<size_t, size_t, size_t>, real> weights_;
        std::map<std::tuple<size_t, size_t, size_t>, gmx::RVec> normals_;



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
                        MolecularPath &molPath);


    private:

        //
        const real PI_ = std::acos(-1.0);

        // 
        inline int numPlanarVertices(real &d, real &r);

        // 
        std::vector<gmx::RVec> generateNormals(
                const std::vector<gmx::RVec> &tangents);
        RegularVertexGrid generateGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                std::vector<SplineCurve1D> &properties,
                std::pair<size_t, size_t> resolution,
                std::pair<real, real> range,
                gmx::RVec chanDirVec);
        void generatePropertyGrid(
                SplineCurve3D &centreLine,
                SplineCurve1D &radius,
                SplineCurve1D &property,
                size_t p,
                gmx::RVec chanDirVec,
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

