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


#ifndef VERSION_HPP
#define VERSION_HPP

#include <string>
#include <iostream>

extern const char g_GIT_SHA1[];
extern const char g_CHAP_VERSION_MAJOR[];
extern const char g_CHAP_VERSION_MINOR[];
extern const char g_CHAP_VERSION_PATCH[];


/*!
 * \brief Returns CHAP version string.
 *
 * The string consists of four point-separated substring representing the 
 * major version number, minor version number, patch number and git hash 
 * respectively.
 */
extern "C" inline std::string
chapVersionString()
{
    return std::string(g_CHAP_VERSION_MAJOR) + std::string(".") + 
           std::string(g_CHAP_VERSION_MINOR) + std::string(".") + 
           std::string(g_CHAP_VERSION_PATCH); 
};


/*!
 * \brief Returns major version number of CHAP.
 */
extern "C" inline std::string
chapVersionMajor()
{
    return std::string(g_CHAP_VERSION_MAJOR);
};


/*!
 * \brief Returns minor version number of CHAP.
 */
extern "C" inline std::string
chapVersionMinor()
{
    return std::string(g_CHAP_VERSION_MINOR);
};


/*!
 * \brief Returns patch number of current version of CHAP.
 */
extern "C" inline std::string
chapVersionPatch()
{
    return std::string(g_CHAP_VERSION_PATCH);
};


/*!
 * \brief Returns git SHA1 hash for current version of CHAP.
 */
extern "C" inline std::string
chapVersionGitHash()
{
    return std::string(g_GIT_SHA1);
};

#endif

