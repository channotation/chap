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


#include <iostream>

#include "config/front_matter.hpp"
#include "config/version.hpp"



/*!
 * Prints a banner and an opening message.
 */
void
FrontMatter::print()
{
	// banner:
	std::cout<<std::endl;
    std::cout<<"		 CCCCCC  HH     HH    AAA    PPPPPPPP  "<<std::endl
             <<"		CC    CC HH     HH   AA AA   PP     PP "<<std::endl
			 <<"		CC       HH     HH  AA   AA  PP     PP "<<std::endl
			 <<"		CC       HHHHHHHHH AA     AA PPPPPPPP  "<<std::endl
			 <<"		CC       HH     HH AAAAAAAAA PP        "<<std::endl
			 <<"		CC    CC HH     HH AA     AA PP        "<<std::endl
			 <<" 		 CCCCCC  HH     HH AA     AA PP        "<<std::endl
			 <<std::endl;

	// name and version info:
    std::cout<<"             "
             <<"The Channel Annotation Package, "
             <<"version "<<chapVersionString()
             <<std::endl<<std::endl;
}

