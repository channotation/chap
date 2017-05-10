#include <cmath>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h> 

#include "io/molecular_path_obj_exporter.hpp"


/*
 * Constructor.
 */
MolecularPathObjExporter::MolecularPathObjExporter()
{
    
}


/*
 *
 */
void
MolecularPathObjExporter::operator()(char *filename,
                                     MolecularPath molPath)
{
    std::cout<<"begin molpath"<<std::endl;

    real d = 0.1;
    real r = 1;

    real extrapDist = 0.0;

    int nLen = 4;
    int nPhi = 10;


    real deltaLen = molPath.length() / (nLen - 1);
    real deltaPhi = 2.0*PI_/(nPhi - 1);

    // get tanges to centre line:
    // ERROR: something goes wrong on lower level here!
    std::vector<gmx::RVec> tangents = molPath.sampleTangents(nLen, extrapDist);

    std::cout<<"post tangent"<<std::endl;

    // construct normals:
    std::vector<gmx::RVec> normals;
    normals.reserve(tangents.size());
    for(unsigned int i = 0; i < tangents.size(); i++)
    {
//        normals.push_back(tangents[i]);
    }
    

    // calculate vertices:
    std::vector<gmx::RVec> vertices;
    vertices.reserve(nLen*nPhi);
    for(int i = 0; i < nLen; i++)
    {
        for(int j = 0; j < nPhi; j++)
        {
//            vertices.push_back(tangents[0]);
        }
    }


    // construct triangular faces:
    std::vector<std::vector<int>> faces;
    for(int i = 0; i < nLen - 1; i++)
    {
        for(int j = 0; j < nPhi; j++)
        {
            // calculate linear indeces:
            int kbl = i*nPhi + j; 
            int kbr = kbl + 1;
            int ktl = kbl + nPhi;
            int ktr = kbr + nPhi;

            // two triangles for each square:
            std::vector<int> faceA = {1, 2, 3};
            std::vector<int> faceB = {1, 2, 3};

            // add to vector of faces:
            faces.push_back(faceA);
            faces.push_back(faceB);   
        }
    }

    // write this to obj file:
    this -> write(filename, vertices, faces);

    std::cout<<"molpath"<<std::endl;
}


/*
 *
 */
int
MolecularPathObjExporter::numPlanarVertices(real &d, real &r)
{
    return std::max(static_cast<int>(std::ceil(PI_/(2.0*std::acos(1.0 - d*d/(2.0*r*r))))), 4);
}

