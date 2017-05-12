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
                                     MolecularPath &molPath)
{
    std::cout<<"begin molpath"<<std::endl;

    real d = 0.1;
    real r = 1;

    real extrapDist = 0.0;

    int nLen = 10;
    int nPhi = 50;


    real deltaLen = molPath.length() / (nLen - 1);
    real deltaPhi = 2.0*PI_/(nPhi - 1);



    // get tangents and points on centre line to centre line:
    std::vector<gmx::RVec> tangents = molPath.sampleTangents(nLen, extrapDist);
    std::vector<gmx::RVec> centrePoints = molPath.samplePoints(nLen, extrapDist);


    // construct orthogonal vector:
    // TODO: this one needs to be updated for every tangent!
    gmx::RVec normal = orthogonalVector(tangents[0]);
 

    //
    std::vector<gmx::RVec> vertices;
    for(unsigned int i = 0; i < tangents.size(); i++)
    {

        // construct sample points:
        for(unsigned int j = 0; j < nPhi; j++)
        {

            // rotate normal vector:
            gmx::RVec rotNormal = rotateAboutAxis(normal, tangents[i], j*deltaPhi);
            std::cout<<"rotNormal = "<<rotNormal[0]<<"  "<<rotNormal[1]<<" "<<rotNormal[2]<<std::endl;

            // generate vertex:
            gmx::RVec vertex = centrePoints[i];
            vertex[XX] += rotNormal[XX];
            vertex[YY] += rotNormal[YY];
            vertex[ZZ] += rotNormal[ZZ];
            vertices.push_back(vertex);
        }
    }











    // construct triangular faces:
    std::vector<std::vector<int>> faces;
    for(int i = 0; i < nLen - 1; i++)
    {
        for(int j = 0; j < nPhi - 1; j++)
        {
            // calculate linear indeces:
            int kbl = i*nPhi + j; 
            int kbr = kbl + 1;
            int ktl = kbl + nPhi;
            int ktr = kbr + nPhi;

            // two triangles for each square:
            std::vector<int> faceA = {kbl + 1, ktl + 1, ktr + 1};
            std::vector<int> faceB = {kbl + 1, ktr + 1, kbr + 1};

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


/*
 *
 */
gmx::RVec
MolecularPathObjExporter::orthogonalVector(gmx::RVec vec)
{
    // find first nonzero element in vector:
    int idxNonZero = -1;
    for(int i = 0; i < 3; i++)
    {
        if( vec[i] > std::numeric_limits<real>::epsilon() )
        {
            idxNonZero = i;
            break;
        }
    }

    // sanity check:
    if( idxNonZero == -1 )
    {
        std::cerr<<"ERROR: Can not find orthogonal to null vector!"<<std::endl;
        std::abort();
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


/*
 *
 */
gmx::RVec
MolecularPathObjExporter::rotateAboutAxis(gmx::RVec vec, 
                                          gmx::RVec axis,
                                          real angle)
{
    // evaluate trigonometrix functions of rotation angle:
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

