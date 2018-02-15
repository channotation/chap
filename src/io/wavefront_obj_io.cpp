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

#include "io/wavefront_obj_io.hpp"


/*!
 * Constructs a face from a vector of vertex indices.
 */
WavefrontObjFace::WavefrontObjFace(
        const std::vector<int> &vertexIdx,
        std::string mtlName)
    : vertexIdx_(vertexIdx)
    , mtlName_(mtlName)
{
    // check uniqueness of indices:
    std::vector<int> tmp = vertexIdx_;
    std::sort(tmp.begin(), tmp.end());
    if( std::adjacent_find(tmp.begin(), tmp.end()) != tmp.end() )
    {
        throw std::logic_error("Vertex indices must be unique.");
    }
}


/*!
 * Constructs a face from vectors of vertex indices and normal indices. These
 * should be of the same size.
 */
WavefrontObjFace::WavefrontObjFace(
        const std::vector<int> &vertexIdx,
        const std::vector<int> &normalIdx,
        std::string mtlName)
{
    // sanity check:
    if( normalIdx.size() != vertexIdx.size() )
    {
        throw std::logic_error("Number of vertex normals must equal number "
                               "of vertices.");
    }

    // assign to internal objects:
    vertexIdx_ = vertexIdx;
    normalIdx_ = normalIdx;
    mtlName_ = mtlName;
}


/*!
 * Returns number of vertices associated with this face.
 */
int
WavefrontObjFace::numVertices() const
{
    return vertexIdx_.size();
}


/*!
 * Returns the i-th vertex index.
 */
int
WavefrontObjFace::vertexIdx(int i) const
{
    return vertexIdx_[i];
}


/*!
 * Returns the i-th normal index.
 */
int
WavefrontObjFace::normalIdx(int i) const
{
    return normalIdx_[i];
}


/*!
 * Returns a flag indicating whether this face has vertex normals.
 */
bool
WavefrontObjFace::hasNormals() const
{
    if( normalIdx_.size() == vertexIdx_.size() )
    {
        return true;
    }
    else
    {
        return false;
    }
}


/*!
 * Construct a new group with given name. Faces can be added using addFace().
 */
WavefrontObjGroup::WavefrontObjGroup(std::string name)
    : groupname_(name)
{
    
}


/*!
 * Copy constructor for allowing use of STL containers.
 */
WavefrontObjGroup::WavefrontObjGroup(const WavefrontObjGroup &other)
    : groupname_(other.groupname_)
    , faces_(other.faces_)
{
    
}


/*!
 * Adds face to internal list of faces. 
 */
void
WavefrontObjGroup::addFace(WavefrontObjFace face)
{
    faces_.push_back(face);
}


/*!
 * Constructs and empty OBJ object with a given name.
 */
WavefrontObjObject::WavefrontObjObject(std::string name)
    : name_(name)
    , mtllib_("")
{
    
}


/*!
 * Add new vertices to the OBJ object. These are appended to the existing list
 * of vertices and there is no check for redundancy. As only geometric vertex
 * coordinates are given, the weight is assumed to be one.
 */
void
WavefrontObjObject::addVertices(const std::vector<gmx::RVec> &vertices)
{
    for(auto vert : vertices)
    {
        // default weight is 1.0:
        std::pair<gmx::RVec, real> vertex(vert, 1.0);
        vertices_.push_back(vertex);
    }
}


/*!
 * Add new vertices to the OBJ object. These are appended to the existing list
 * of vertices and there is no check for redundancy.  
 */
void
WavefrontObjObject::addVertices(
        const std::vector<std::pair<gmx::RVec, real>> &vertices)
{
    for(auto vert : vertices)
    {
        // sanity check:
        if( vert.second < 0.0 || vert.second > 1.0 )
        {
            throw std::logic_error("Wavefront OBJ vertex weights must be in "
                                   "unit interval.");
        }

        // add to vertex storage:
        vertices_.push_back(vert);
    }
}


/*!
 * Add new vertex normals to the OBJ object. These are appended to the existing
 * list of vertex normals and the user is responsible for ensuring that the 
 * final number of vertex normals is either equal to the number of vertices or
 * zero.
 */
void
WavefrontObjObject::addVertexNormals(const std::vector<gmx::RVec> &normals)
{
    // accept normals:
    normals_.insert(normals_.end(), normals.begin(), normals.end());
}


/*!
 * Adds a new named group to the OBJ object.
 */
void
WavefrontObjObject::addGroup(const WavefrontObjGroup &group)
{
    groups_.push_back( group );
}


/*!
 * Sets the name of the material library to be used. Input argument should 
 * include the .mtl extension.
 */
void
WavefrontObjObject::setMaterialLibrary(std::string mtl)
{
    mtllib_ = mtl;
}


/*!
 * Scales the shape by a given factor. To this end, all vertex positions are 
 * shifted so that the centre of geometry is the origin. Then all position 
 * vectors are multiplied by a given factor before finally the positions are
 * shifted back again with the scaling factor also applied to the the shift
 * vector.
 */
void
WavefrontObjObject::scale(real fac)
{
    // shift vertices to be centred around origin:
    gmx::RVec shift = calculateCog();
    shift[XX] = -shift[XX];
    shift[YY] = -shift[YY];
    shift[ZZ] = -shift[ZZ];
    this -> shift(shift);

    // scale all position vectors:
    for(auto it = vertices_.begin(); it != vertices_.end(); it++)
    {
        (*it).first[XX] *= fac;
        (*it).first[YY] *= fac;
        (*it).first[ZZ] *= fac;
    }

    // shift vertices back to original centre of geometry:
    shift[XX] *= -fac;
    shift[YY] *= -fac;
    shift[ZZ] *= -fac;
    this -> shift(shift);
}


/*!
 * Shifts all vertex positions by the given vector.
 */
void WavefrontObjObject::shift(gmx::RVec shift)
{
    for(auto it = vertices_.begin(); it != vertices_.end(); it++)
    {
        (*it).first[XX] += shift[XX];
        (*it).first[YY] += shift[YY];
        (*it).first[ZZ] += shift[ZZ];
    }    
}


/*!
 * Calculates the centre of geometry of all vertices.
 */
gmx::RVec
WavefrontObjObject::calculateCog()
{
    gmx::RVec cog(0.0, 0.0, 0.0);

    for(auto it = vertices_.begin(); it != vertices_.end(); it++)
    {
        cog[XX] += (*it).first[XX];
        cog[YY] += (*it).first[YY];
        cog[ZZ] += (*it).first[ZZ];
    }

    cog[XX] /= vertices_.size();
    cog[YY] /= vertices_.size();
    cog[ZZ] /= vertices_.size();

    return cog;
}


/*!
 * Returns a flag indicating if the Wavefront object is valid. In particular,
 * it check that each vertex or vertex normal referenced by the faces is 
 * present in the vertex or vertex normal set.
 */
bool
WavefrontObjObject::valid() const
{
    // loop over groups:
    for( auto& group : groups_ )
    {
        // loop over faces:
        for( auto& face : group.faces_ )
        {
            // loop over vertices:
            for( size_t i = 0; i < face.numVertices(); i++ )
            {
                // are all referenced vertices present?
                if( face.vertexIdx(i) > vertices_.size() ||
                    face.vertexIdx(i) < 1 )
                {
                    return false;
                }

                // are all referenced normals present?
                if( face.hasNormals() )
                {
                    if( face.normalIdx(i) > normals_.size() ||
                        face.normalIdx(i) < 1 )
                    {
                        return false;
                    }
                }
            }
        }
    }

    // if nothing failed, return true:
    return true;
}


/*!
 * Writes an OBJ object to a file of the given name. 
 */
void
WavefrontObjExporter::write(std::string fileName,
                            WavefrontObjObject object)
{
    // sanity checks:
    if( !object.valid() )
    {
        throw std::logic_error("WavefrontObjExporter encountered invalid "
                               "OBJ object.");
    }

    // open file stream:
    obj_.open(fileName.c_str(), std::fstream::out);

    // writer header comment:
    writeComment("produced by CHAP");

    // write name of material library:
    writeMaterialLibrary(object.mtllib_);

    // write object name:
    writeObject(object.name_);

    // write vertices:
    obj_<<std::endl;
    for(unsigned int i = 0; i < object.vertices_.size(); i++)
    {
        writeVertex(object.vertices_[i]);
    }

    // write vertex normals:
    obj_<<std::endl;
    for(unsigned int i = 0; i < object.normals_.size(); i++)
    {
        writeVertexNormal(object.normals_[i]);
    }

    // write groups:
    std::vector<WavefrontObjGroup>::iterator it;
    for(it = object.groups_.begin(); it != object.groups_.end(); it++)
    {
        // write group name:
        writeGroup(it -> groupname_);

        // write faces in this group:
        for(size_t j = 0; j < it -> faces_.size(); j++)
        {
            writeFace(it -> faces_[j]);
        }
    }

    // close file stream:
    obj_.close();
}


/*!
 * Writes a comment line to an OBJ file.
 */
void
WavefrontObjExporter::writeComment(std::string comment)
{
    obj_ <<"# "<<comment<<std::endl;
}


/*!
 * Writes material library referenct to an OBJ file.
 */
void
WavefrontObjExporter::writeMaterialLibrary(std::string mtl)
{
    // library set?
    if( mtl != "" )
    {
        obj_<<"mtllib "<<mtl<<std::endl;
    }
}


/*!
 * Writes a group line to an OBJ file.
 */
void
WavefrontObjExporter::writeGroup(std::string group)
{
    obj_ <<std::endl;
    obj_ <<"g "<<group<<std::endl;
}


/*!
 * Writes an object line to an OBJ file.
 */
void
WavefrontObjExporter::writeObject(std::string object)
{
    obj_ <<std::endl;
    obj_ <<"o "<<object<<std::endl;
}


/*!
 * Writes a vertex entry to an OBJ file.
 */
void
WavefrontObjExporter::writeVertex(std::pair<gmx::RVec, real> vertex)
{
    obj_<<"v "<<vertex.first[XX]<<" "
              <<vertex.first[YY]<<" "
              <<vertex.first[ZZ]<<" "
              <<vertex.second<<std::endl;
}


/*!
 * Writes a vertex normal to an OBJ file.
 */
void
WavefrontObjExporter::writeVertexNormal(gmx::RVec norm)
{
    obj_<<"vn "<<norm[XX]<<" "
               <<norm[YY]<<" "
               <<norm[ZZ]<<std::endl;
}


/*
 *
 */
void
WavefrontObjExporter::writeFace(const WavefrontObjFace &face)
{
    // need to write material?
    if( face.mtlName_ != "" )
    {
        // different from current material?
        if( face.mtlName_ != crntMtlName_ )
        {
            crntMtlName_ = face.mtlName_;
            obj_<<"usemtl "<<face.mtlName_<<std::endl;
        }
    }

    // write actual face entry:
    obj_<<"f ";

    for(size_t i = 0; i < face.numVertices(); i++)
    {
        obj_<<face.vertexIdx(i);

        if( face.hasNormals() )
        {
            obj_<<"//"<<face.normalIdx(i);
        }

        obj_<<" ";
    }

    obj_<<std::endl;
}

