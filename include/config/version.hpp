#ifndef VERSION_HPP
#define VERSION_HPP

#include <string>
#include <iostream>

const int test = 1;

extern const char g_GIT_SHA1[];
extern const char g_CHAP_VERSION_MAJOR[];
extern const char g_CHAP_VERSION_MINOR[];
extern const char g_CHAP_VERSION_PATCH[];



/*
 *
 */
extern "C" inline std::string
chapVersionString()
{
    return std::string(g_CHAP_VERSION_MAJOR) + std::string(".") + 
           std::string(g_CHAP_VERSION_MINOR) + std::string(".") + 
           std::string(g_CHAP_VERSION_PATCH) + std::string(".") + 
           std::string(g_GIT_SHA1);
};

#endif

