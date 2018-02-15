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


#ifndef WAVEFRONT_MTL_IO_HPP
#define WAVEFRONT_MTL_IO_HPP

#include <fstream>
#include <string>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h>


/*!
 * \brief Representation of material in MTL format.
 *
 * This is an abstract container for material specification in the Wavefront
 * MTL format. 
 *
 * \note Currently only a small subset of material properties is supported.
 */
class WavefrontMtlMaterial
{
    friend class WavefrontMtlExporter;

    public:

        // constructor::
        WavefrontMtlMaterial(std::string name);

        // setter functions for material properties:
        void setAmbientColour(gmx::RVec col);
        void setDiffuseColour(gmx::RVec col);
        void setSpecularColour(gmx::RVec col);

    private:

        // material properties:
        std::string name_;
        gmx::RVec ambientColour_;
        gmx::RVec diffuseColour_;
        gmx::RVec specularColour_;

};


/*!
 * \brief Abstract representation of an MTL file.
 *
 * This is essentially a container for WavefrontMtlMaterials that serves to
 * unify the input to WavefrontMtlExporter.
 */
class WavefrontMtlObject
{
    friend class WavefrontMtlExporter;

    public:

        // interface for adding materials:
        void addMaterial(WavefrontMtlMaterial material);

    private:

        // collection of materials:
        std::vector<WavefrontMtlMaterial> materials_;
};


/*!
 * \brief Serialiser for material specification in Wavefront MTL format.
 *
 * This class serves to write the material specifications contained in a
 * WavefrontMtlObject to file.
 */
class WavefrontMtlExporter
{
    public:

        // public interfacr for writing material definitions to file:
        void write(std::string fileName, WavefrontMtlObject object);

    private:

        // file handle:
        std::fstream file_;
    
        // utilities for writing individual lines:
        void writeMaterialName(std::string name);
        void writeAmbientColour(const gmx::RVec &col);
        void writeDiffuseColour(const gmx::RVec &col);
        void writeSpecularColour(const gmx::RVec &col);
};

#endif

