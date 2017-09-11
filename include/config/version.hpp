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

