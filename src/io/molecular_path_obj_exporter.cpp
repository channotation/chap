#include <cmath>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h> 

#include "io/molecular_path_obj_exporter.hpp"
#include "io/wavefront_obj_io.hpp"


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
MolecularPathObjExporter::operator()(std::string fileName,
                                     MolecularPath &molPath)
{

//    real d = 0.1;
//    real r = 1;

    real extrapDist = 1.0;

    int nLen = 100;
    int nPhi = 100;


//    real deltaLen = molPath.length() / (nLen - 1);
    real deltaPhi = 2.0*PI_/(nPhi - 1);



    // get tangents and points on centre line to centre line:
    std::vector<gmx::RVec> tangents = molPath.sampleNormTangents(nLen, extrapDist);

    std::vector<gmx::RVec> centrePoints = molPath.samplePoints(nLen, extrapDist);
    std::vector<real> radii = molPath.sampleRadii(nLen, extrapDist);

    // preallocate output vertices:
    std::vector<gmx::RVec> vertices;
    vertices.reserve(tangents.size()*nPhi);


    // construct vertices around first point:
    gmx::RVec normal = orthogonalVector(tangents[0]);
    unitv(normal, normal);
    std::vector<gmx::RVec> newVertices = vertexRing(centrePoints[0],
                                                    tangents[0],
                                                    normal,
                                                    radii[0],
                                                    deltaPhi,
                                                    nPhi);
    vertices.insert(vertices.end(), newVertices.begin(), newVertices.end());


    // follow spline curve and construct vertices around all other points:
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

        // construct sample points:
        newVertices = vertexRing(centrePoints[i],
                                 tangents[i],
                                 normal,
                                 radii[i],
                                 deltaPhi,
                                 nPhi);
        vertices.insert(vertices.end(), newVertices.begin(), newVertices.end());
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

    // create an OBJ object:
    WavefrontObjObject obj("pore");
    obj.addVertices(vertices);
    obj.addGroup("inner_surface", faces);

    // scale object by factor of 10 to convert nm to Ang:
    obj.scale(10.0);
    obj.calculateCog();

    // create OBJ exporter and write to file:
    WavefrontObjExporter objExp;
    objExp.write(fileName, obj);
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
std::vector<gmx::RVec> 
MolecularPathObjExporter::vertexRing(gmx::RVec base,
                                     gmx::RVec tangent,
                                     gmx::RVec normal,
                                     real radius,
                                     real angleIncrement,
                                     size_t nIncrements)
{
    // make sure input normal vector is normalised:
    unitv(normal, normal);

    // preallocate vertex vector:
    std::vector<gmx::RVec> vertices;
    vertices.reserve(nIncrements);

    // sample vertices in a ring around the base point: 
    for(size_t j = 0; j < nIncrements; j++)
    {
        // rotate normal vector:
        gmx::RVec rotNormal = rotateAboutAxis(normal, tangent, j*angleIncrement);

        // generate vertex:
        gmx::RVec vertex = base;
        vertex[XX] += radius*rotNormal[XX];
        vertex[YY] += radius*rotNormal[YY];
        vertex[ZZ] += radius*rotNormal[ZZ];
        vertices.push_back(vertex);
    }

    // return ring of vertices:
    return vertices;
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
        if( std::abs(vec[i]) > std::numeric_limits<real>::epsilon() )
        {
            idxNonZero = i;
            break;
        }
    }

    // sanity check:
    if( idxNonZero == -1 )
    {
        std::cerr<<"ERROR: Can not find orthogonal to null vector!"<<std::endl;
        std::cerr<<"vec = "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<std::endl;
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

