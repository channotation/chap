#include <algorithm>
#include <cmath>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h> 

#include "io/molecular_path_obj_exporter.hpp"


/*
 *
 */
RegularVertexGrid::RegularVertexGrid(
        std::vector<real> s,
        std::vector<real> phi)
    : s_(s)
    , phi_(phi)
{
    
}


/*
 *
 */
void
RegularVertexGrid::addVertex(size_t i, size_t j, gmx::RVec vertex, real weight)
{
    std::pair<size_t, size_t> key(i, j);
    vertices_[key] = vertex;
    weights_[key] = weight;
}


/*
 *
 */
void
RegularVertexGrid::addVertexNormal(size_t i, size_t j, gmx::RVec normal)
{
    std::pair<size_t, size_t> key(i, j);
    normals_[key] = normal;
}







/*
 *
 */
std::vector<gmx::RVec>
RegularVertexGrid::vertices()
{
    // build linearly indexed vector from dual indexed map:
    std::vector<gmx::RVec> vert;
    vert.reserve(vertices_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::pair<size_t, size_t> key(i, j);
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


/*
 *
 */
std::vector<std::pair<gmx::RVec, real>>
RegularVertexGrid::weightedVertices()
{
    // build linearly indexed vector from dual indexed map:
    std::vector<std::pair<gmx::RVec, real>> vert;
    vert.reserve(vertices_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::pair<size_t, size_t> key(i, j);
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


/*
 *
 */
void
RegularVertexGrid::normalsFromFaces()
{
    // array sizes for index wrap:
    int mI = s_.size();
    int mJ = phi_.size();

    for(int i = 0; i < s_.size(); i++)
    {
        for(int j = 0; j < phi_.size(); j++)
        {
            // index pairs for neighbouring vertices:
            // TODO: not circular in i, cant wrap around
            std::pair<size_t, size_t> crntKey(i, j);
            std::pair<size_t, size_t> leftKey(i, (j - 1 + mJ) % mJ);
            std::pair<size_t, size_t> rghtKey(i, (j + 1 + mJ) % mJ);
            std::pair<size_t, size_t> upprKey((i + 1 + mI) % mI, j);
            std::pair<size_t, size_t> lowrKey((i - 1 + mI) % mI, j);
            std::pair<size_t, size_t> dglrKey((i - 1 + mI) % mI, (j + 1 + mJ) % mJ);
            std::pair<size_t, size_t> dgulKey((i + 1 + mI) % mI, (j - 1 + mJ) % mJ);

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


/*
 *
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


/*
 *
 */
std::vector<gmx::RVec>
RegularVertexGrid::normals()
{
    std::vector<gmx::RVec> norm;
    norm.reserve(normals_.size());
    for(size_t i = 0; i < s_.size(); i++)
    {
        for(size_t j = 0; j < phi_.size(); j++)
        {
            std::pair<size_t, size_t> key(i, j);
            if( normals_.find(key) != normals_.end() )
            {
                norm.push_back(normals_[key]); 
            }
            else
            {
//                throw std::logic_error("Invalid vertex normal reference "
//                                       "encountered.");
            }
        }
    }
    return norm;
}





/*
 *
 */
std::vector<WavefrontObjFace>
RegularVertexGrid::faces()
{
    // sanity checks:
    if( phi_.size() * s_.size() != vertices_.size() )
    {
        std::cout<<"phi.size = "<<phi_.size()<<"  "
                 <<"s.size = "<<s_.size()<<"  "
                 <<"vert.size = "<<vertices_.size()<<"  "
                 <<std::endl;
        throw std::logic_error("Cannot generate faces on incomplete grid.");
    }

    if( !normals_.empty() && normals_.size() != vertices_.size() )
    {
        throw std::logic_error("Number of vertex normals does not equal "
                               "number of vertices.");
    }

    // preallocate face vector:
    std::vector<WavefrontObjFace> faces;
    faces.reserve(phi_.size()*s_.size());

    // loop over grid:
    for(size_t i = 0; i < s_.size() - 1; i++)
    {
        for(size_t j = 0; j < phi_.size() - 1; j++)
        {
            // calculate linear indices:
            int kbl = i*phi_.size() + j; 
            int kbr = kbl + 1;
            int ktl = kbl + phi_.size();
            int ktr = kbr + phi_.size();

            // two faces per square:
            if( normals_.empty() )
            {
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, ktr + 1, ktl + 1}) );
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, kbr + 1, ktr + 1}) );
            }
            else
            {
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, ktr + 1, ktl + 1},
                        {kbl + 1, ktr + 1, ktl + 1}) );
                faces.push_back( WavefrontObjFace(
                        {kbl + 1, kbr + 1, ktr + 1},
                        {kbl + 1, kbr + 1, ktr + 1}) );
            }
        }
    }

    // wrap around:
    for(size_t i = 0; i < s_.size() - 1; i++)
    {
        // calculate linear indices:
        int kbl = i*phi_.size() + phi_.size() - 1; 
        int kbr = i*phi_.size();
        int ktl = kbl + phi_.size();
        int ktr = kbr + phi_.size();

        // two faces per square:
        if( normals_.empty() )
        {
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, ktr + 1, ktl + 1}) );
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, kbr + 1, ktr + 1}) );
        }
        else
        {

            faces.push_back( WavefrontObjFace(
                    {kbl + 1, ktr + 1, ktl + 1},
                    {kbl + 1, ktr + 1, ktl + 1}) );
            faces.push_back( WavefrontObjFace(
                    {kbl + 1, kbr + 1, ktr + 1},
                    {kbl + 1, kbr + 1, ktr + 1}) );
        }
    }

    // return face vector:
    return faces;
}



/*
 * Constructor.
 */
MolecularPathObjExporter::MolecularPathObjExporter()
{
    
}


void
MolecularPathObjExporter::operator()(std::string fileName,
                                     MolecularPath &molPath)
{
    // define evaluation range:   
    real extrapDist = 0.0;
    std::pair<real, real> range(molPath.sLo() - extrapDist,
                                molPath.sHi() + extrapDist);

    // define resolution:
    int numPhi = 30;
    int numLen = std::pow(2, 8) + 1;
    std::pair<size_t, size_t> resolution(numLen, numPhi);
    
    // spline curves for variuous properties:
    auto centreLine = molPath.centreLine();
    auto pathRadius = molPath.pathRadius();

    
    gmx::RVec chanDirVec(0.0, 0.0, 1.0);
    RegularVertexGrid grid = generateGrid(
            centreLine,
            pathRadius,
            pathRadius,
            resolution,
            range,
            chanDirVec);

    // obtain vertices, normals, and faces from grid:
    grid.normalsFromFaces();
    auto vertices = grid.weightedVertices();
    auto vertexNormals = grid.normals();
    auto faces = grid.faces();




    std::cout<<"numVert = "<<vertices.size()<<"  "
             <<"normNorm = "<<vertexNormals.size()<<"  "
             <<std::endl;

    

    // add faces to surface:
    WavefrontObjGroup surface("pore_surface");
    for(auto face : faces)
    {
/*        std::cout<<"face = "<<"  "
                 <<face.vertexIdx_[0]<<"  "
                 <<face.vertexIdx_[1]<<"  "
                 <<face.vertexIdx_[2]<<"  "
                 <<"hasNormals = "<<face.hasNormals()<<"  "
                 <<std::endl;*/
        surface.addFace(face);
    }
   

    // assemble OBJ object:
    WavefrontObjObject obj("pore");
    obj.addVertices(vertices);
    obj.addVertexNormals(vertexNormals);
    obj.addGroup(surface);

    // scale object by factor of 10 to convert nm to Ang:
    obj.scale(10.0);
    obj.calculateCog();

    // create OBJ exporter and write to file:
    WavefrontObjExporter objExp;
    objExp.write(fileName, obj);

    
}


/*!
 * Generates grid along MolecularPath, but uses a tangent vector equal to the
 * channel direction vector. This does not strictly give the correct surface, 
 * as the tangent may vary along the centre line. However, this method will not
 * create overlapping vertices, so long as the centre line does not change 
 * direction with respect to the channel direction vectors, which is guaranteed
 * for the inplane optimised probe path finder. Moreover, all surface points 
 * generated this way are guaranteed to not lie outside the pore, so long as
 * a probe based method was used (this is true because the free distance around 
 * a centre is free in any direction). 
 */
RegularVertexGrid
MolecularPathObjExporter::generateGrid(
        SplineCurve3D &centreLine,
        SplineCurve1D &radius,
        SplineCurve1D &property,
        std::pair<size_t, size_t> resolution,
        std::pair<real, real> range,
        gmx::RVec chanDirVec)
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

    // normalise channel direction vector:
    unitv(chanDirVec, chanDirVec);

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

    // sample points, radii, and tangents along molecular path:
    std::vector<gmx::RVec> centres;
    centres.reserve(numLen);
    std::vector<real> radii;
    radii.reserve(numLen);
    std::vector<real> prop;
    prop.reserve(numLen);
    for(auto eval : s)
    {
        centres.push_back( centreLine.evaluate(eval, 0) );
        radii.push_back( radius.evaluate(eval, 0) );
        prop.push_back( property.evaluate(eval, 0) );
    }

    // sample scalar property along the path:
    shiftAndScale(prop);

    // tangents always in direction of channel:
    std::vector<gmx::RVec> tangents(centres.size(), chanDirVec);

    // sample normals along molecular path:
    auto normals = generateNormals(tangents);
    
    // store vertex rings:
    std::vector<gmx::RVec> vertRing(phi.size());

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

        // add to grid:
        grid.addVertex(idxLen, k, vertex, prop[idxLen]);
    }

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

        // add to grid:
        grid.addVertex(idxLen, k, vertex, prop[idxLen]);
    }

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

                vertRing[k] = vertex;
            }
 
            // add vertices to grid:
            for(size_t k = 0; k < vertRing.size(); k++)
            {
                grid.addVertex(idxLen, k, vertRing[k], prop[idxLen]);
            }
        }
    }

    // return grid:
    return grid;
}


/*
 *
 */
RegularVertexGrid
MolecularPathObjExporter::generateGrid(
        MolecularPath &molPath,
        size_t numLen,
        size_t numPhi,
        real extrapDist)
{    
    // check if number of intervals is power of two:
    int numInt = numLen - 1;
    if( (numInt & (numInt)) == 0 && numInt > 0 )
    {
        throw std::logic_error("Number of steps along pore must be power of "
                               "two.");
    }

    // generate grid coordinates:
    auto s = molPath.sampleArcLength(numLen, extrapDist);
    std::vector<real> phi;
    phi.reserve(numPhi);
    for(size_t i = 0; i < numPhi; i++)
    {
        phi.push_back(i*2.0*M_PI/numPhi);
    }

    // generate grid from coordinates:
    RegularVertexGrid grid(s, phi);

    // sample centre points radii, and tangents along molecular path:
    auto centres = molPath.samplePoints(s);
//    auto radii  = molPath.sampleRadii(s);
//    auto tangents = molPath.sampleNormTangents(s);
    std::vector<gmx::RVec> tangents(centres.size(), gmx::RVec(0.0, 0.0, 1.0));

    std::vector<real> radii(tangents.size(), 0.025);
    std::vector<real> prop(tangents.size(), 1.0);
    
    // sample normals along molecular path:
    auto normals = generateNormals(tangents);
    

    // store vertex rings:
    std::vector<gmx::RVec> vertRing(phi.size());
    std::vector<std::vector<gmx::RVec>> vertexRings(tangents.size());


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

        // add to grid:
        grid.addVertex(idxLen, k, vertex, prop[idxLen]);
        grid.addVertexNormal(idxLen, k, rotNormal);
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

        // add to grid:
        grid.addVertex(idxLen, k, vertex, prop[idxLen]);
        grid.addVertexNormal(idxLen, k, rotNormal);
    }
    vertexRings[idxLen] = vertRing;

    // build intermediate vertex rings:
    int numClash = 0;
    int clashLower = 0;
    int clashUpper = 0;
    int numVert = 0;
    int numSpecial = 0;
    int numPers = 0;
    for(int i = 1; i <= numLen; i *= 2)
    {
        for(int j = 1; j < i; j += 2)
        {
            idxLen = j*(numLen - 1)/i;
            int idxLower = (j - 1)*(numLen - 1)/i;
            int idxUpper = (j + 1)*(numLen - 1)/i;

            // FIXME
            tangents[idxLen] = gmx::RVec(0, 0, 1);
            tangents[idxLower] = gmx::RVec(0, 0, 1);
            tangents[idxUpper] = gmx::RVec(0, 0, 1);

            // create ring of vertices:
            bool hasClashes = false;
            std::vector<gmx::RVec> vertNormRing;
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

                // distance vector from upper and lower circle centres:
                gmx::RVec distLower;
                rvec_sub(vertex, centres[idxLower], distLower);
                gmx::RVec distUpper;
                rvec_sub(vertex, centres[idxUpper], distUpper);

                // project distance vector onto tangent:
                real cosLower = iprod(tangents[idxLower], distLower);                
                real cosUpper = iprod(tangents[idxUpper], distUpper);                


                real angleThreshold = 0.5;
                if( cosUpper > angleThreshold || cosLower < angleThreshold )
                {
                    hasClashes = true;
                }

                vertRing[k] = vertex;


                // BELOW: second attempt to fix
                ///////////////////////////////////////////////////////////////


/*
                // distance vector from upper and lower circle centres:
                gmx::RVec distLower;
                rvec_sub(vertex, centres[idxLower], distLower);
                gmx::RVec distUpper;
                rvec_sub(vertex, centres[idxUpper], distUpper);

                // project distance vector onto tangent:
                real cosLower = iprod(tangents[idxLower], distLower);                
                real cosUpper = iprod(tangents[idxUpper], distUpper);                
                numVert++;

                real angleThreshold = 0.0;


                if( cosLower < angleThreshold && cosUpper > angleThreshold )
                {
                    numSpecial++;
                    std::cout<<"SPECIAL CASE"<<std::endl;
//                    std::abort();
                }

                real test = iprod(tangents[idxLower], tangents[idxUpper]);

 
                // check for clashes:
                if( cosLower < angleThreshold )
                {
               std::cout<<"LOWER clash"
                        <<"numRings = "<<vertexRings.size()<<"  "
                         <<"idxLen = "<<idxLen<<"  "
                         <<"idxLower = "<<idxLower<<"  "
                         <<"idxUpper = "<<idxUpper<<"  "
//                         <<"lower.size = "<<vertexRings[idxLower].size()<<"  "
//                         <<"upper.size = "<<vertexRings[idxUpper].size()<<"  "
                         <<"cosLower = "<<cosLower<<"  "
                         <<"cosUpper = "<<cosUpper<<"  "
                         <<"dLo = "<<norm(distLower)<<"  "
                         <<"dUp = "<<norm(distUpper)<<"  "
                         <<"test = "<<test<<"  "
//                         <<"cux = "<<centres[idxUpper][XX]<<"  "
//                         <<"cuy = "<<centres[idxUpper][YY]<<"  "
//                         <<"cuz = "<<centres[idxUpper][ZZ]<<"  "
//                         <<"vx = "<<vertex[XX]<<"  "
//                         <<"vy = "<<vertex[YY]<<"  "
//                         <<"vz = "<<vertex[ZZ]<<"  "
                         <<std::endl;
                    clashLower++;
                    // reuse lower vertex:
                    gmx::RVec vert;
                    rvec_add(
                            vertexRings.at(idxLower).at(k), 
                            vertexRings.at(idxUpper).at(k), 
                            vert);
                    svmul(0.5, vert, vert);

//                    vertRing[k] = vertexRings.at(idxLower).at(k);
                    vertRing[k] = vertex;
                }
                else if( cosUpper > angleThreshold )
                {
               std::cout<<"UPPER clash"
                        <<"numRings = "<<vertexRings.size()<<"  "
                         <<"idxLen = "<<idxLen<<"  "
                         <<"idxLower = "<<idxLower<<"  "
                         <<"idxUpper = "<<idxUpper<<"  "
//                         <<"lower.size = "<<vertexRings[idxLower].size()<<"  "
//                         <<"upper.size = "<<vertexRings[idxUpper].size()<<"  "
                         <<"cosLower = "<<cosLower<<"  "
                         <<"cosUpper = "<<cosUpper<<"  "
                         <<"dLo = "<<norm(distLower)<<"  "
                         <<"dUp = "<<norm(distUpper)<<"  "
                         <<"test = "<<test<<"  "
//                         <<"cux = "<<centres[idxUpper][XX]<<"  "
//                         <<"cuy = "<<centres[idxUpper][YY]<<"  "
//                         <<"cuz = "<<centres[idxUpper][ZZ]<<"  "
//                         <<"vx = "<<vertex[XX]<<"  "
//                         <<"vy = "<<vertex[YY]<<"  "
//                         <<"vz = "<<vertex[ZZ]<<"  "
                         <<std::endl;
                    clashUpper++;
                    gmx::RVec vert;
                    rvec_add(
                            vertexRings.at(idxLower).at(k), 
                            vertexRings.at(idxUpper).at(k), 
                            vert);
                    svmul(0.5, vert, vert);
                    // reuse upper vertex:
//                    vertRing[k] = vertexRings.at(idxUpper).at(k);
                    vertRing[k] = vertex;
                     
                }
                else
                {
                    // use new vertex:
                    vertRing[k] = vertex;
                }
*/




                vertNormRing.push_back(rotNormal);
            }
 

/*
            if( hasClashes == true )
            {
                std::cout<<"has clashes!"<<std::endl;

                // form averaged tangent:
                rvec_add(tangents[idxLower], tangents[idxUpper], tangents[idxLen]);
                svmul(0.5, tangents[idxLen], tangents[idxLen]);
                unitv(tangents[idxLen], tangents[idxLen]);

                // recompute vertices:
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

                    vertRing[k] = vertex;
                }
            }
*/














            for(auto vert : vertRing)
            {
                // distance vector from upper and lower circle centres:
                gmx::RVec distLower;
                rvec_sub(vert, centres[idxLower], distLower);
                gmx::RVec distUpper;
                rvec_sub(vert, centres[idxUpper], distUpper);

                // project distance vector onto tangent:
                real cosLower = iprod(tangents[idxLower], distLower);                
                real cosUpper = iprod(tangents[idxUpper], distUpper);                

                real angleThreshold = 0.0;
                if( cosLower < angleThreshold )
                {
                    std::cout<<"PERSISTENT LOWER CLASH"<<std::endl;
                    numPers++;
                }
                if( cosUpper > angleThreshold )
                {
                    std::cout<<"PERSISTENT Upper CLASH"<<std::endl;
                    numPers++;
                }
            }


            


            // add vertices to grid:
            for(size_t k = 0; k < vertRing.size(); k++)
            {
                grid.addVertex(idxLen, k, vertRing[k], prop[idxLen]);
            }
 
            // add vertex normals to grid:
            for(size_t k = 0; k < vertNormRing.size(); k++)
            {
                grid.addVertexNormal(idxLen, k, vertNormRing[k]);
            }

            // add vertices to ring storage:
            vertexRings[idxLen] = vertRing;


            // store 


            std::cout<<"i = "<<i<<"  "
                     <<"j = "<<j<<"  "
                     <<"numLen = "<<numLen<<"  "
                     <<"crnt = "<<idxLen<<"  "
                     <<"s.size = "<<s.size()<<"  "
                     <<"lo = "<<(j-1)*(numLen - 1)/i<<"  "
                     <<"hi = "<<(j+1)*(numLen - 1)/i<<"  "
                     <<std::endl;
        }
    }

    std::cout<<"numClash = "<<numClash<<std::endl;


    std::cout<<"clashUpper = "<<clashUpper<<"  "
             <<"clashLower = "<<clashLower<<"  "
             <<"numVert = "<<numVert<<"  "
             <<"numSpecial = "<<numSpecial<<"  "
             <<"numPers = "<<numPers<<"  "
             <<std::endl;

    // correct mesh:
    
    // FIXME I think the porblem is that interpolation in XYZ does also shift
    // the vertex in S direction. Might be better to interpolate R if a con-
    // flict is found...
//    grid.interpolateMissing();

    // return grid:
    return grid;
}



/*
 *
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
std::pair<std::vector<gmx::RVec>, std::vector<gmx::RVec>>
MolecularPathObjExporter::vertexRing(gmx::RVec base,
                                     gmx::RVec tangent,
                                     gmx::RVec normal,
                                     real radius,
                                     real angleIncrement,
                                     size_t nIncrements)
{
    // make sure input normal vector is normalised:
    unitv(normal, normal);

    // preallocate vertex and normal vectors:
    std::vector<gmx::RVec> vertices;
    vertices.reserve(nIncrements);
    std::vector<gmx::RVec> normals;
    vertices.reserve(nIncrements);

    // sample vertices in a ring around the base point: 
    for(size_t j = 0; j < nIncrements; j++)
    {
        // rotate normal vector:
        gmx::RVec rotNormal = rotateAboutAxis(normal, tangent, j*angleIncrement);
        normals.push_back(rotNormal);

        // generate vertex:
        gmx::RVec vertex = base;
        vertex[XX] += radius*rotNormal[XX];
        vertex[YY] += radius*rotNormal[YY];
        vertex[ZZ] += radius*rotNormal[ZZ];
            std::cout<<"vex = "<<vertex[XX]<<"  "
                     <<"vey = "<<vertex[YY]<<"  "
                     <<"vez = "<<vertex[ZZ]<<"  "
                     <<"bx = "<<base[XX]<<"  "
                     <<"by = "<<base[YY]<<"  "
                     <<"bz = "<<base[ZZ]<<"  "
                     <<"rx = "<<rotNormal[XX]<<"  "
                     <<"ry = "<<rotNormal[YY]<<"  "
                     <<"rz = "<<rotNormal[ZZ]<<"  "
                     <<"radius = "<<radius<<"  "
                     <<std::endl;
        vertices.push_back(vertex);
    }

    // return ring of vertices:
    std::pair<std::vector<gmx::RVec>, std::vector<gmx::RVec>> ret;
    ret.first = vertices;
    ret.second = normals;
    return ret;
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


/*
 *
 *
 */
real
MolecularPathObjExporter::cosAngle(
        const gmx::RVec &vecA,
        const gmx::RVec &vecB)
{
    return iprod(vecA, vecB) / ( norm(vecA) * norm(vecB) );
}


/*!
 * Shift ands scales all values in input vector so that they lie in the unit
 * interval.
 */
void
MolecularPathObjExporter::shiftAndScale(std::vector<real> &prop)
{
    // find data range:
    real minProp = *std::min_element(prop.begin(), prop.end()); 
    real maxProp = *std::max_element(prop.begin(), prop.end());

    // mind special case of constant property value:
    real shift = -minProp + 1.0;
    real scale = 1.0;
    if( std::fabs(maxProp - minProp) > std::numeric_limits<real>::epsilon() )
    {
        // shift and scaling factor:
        shift = -minProp;
        scale = 1.0/(maxProp - minProp); 
    }

    // shift and scale the property array:
    std::for_each(prop.begin(), prop.end(), [shift](real &p){p += shift;});
    std::for_each(prop.begin(), prop.end(), [scale](real &p){p *= scale;});
}

