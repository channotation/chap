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

#include "config/back_matter.hpp"


/*!
 * Prints citation message.
 */
void
BackMatter::print()
{
    std::cout<<std::endl<<std::endl;
    std::cout<<"Thank you for using the Channel Annotation Package! "
             <<"Please cite:"<<std::endl<<std::endl
             <<"Jemma L. Trick, Sivapalan Chelvaniththilan, Gianni Klesse, Prafulla Aryal, "<<std::endl
             <<"E. Jayne Wallace, Stephen J. Tucker,Mark S.P. Sansom"<<std::endl
             <<"Functional Annotation of Ion Channel Structures by Molecular Simulation"<<std::endl
             <<"In Structure, Volume 24, Issue 12, 2016, Pages 2207-2216, ISSN 0969-2126"<<std::endl
             <<"https://doi.org/10.1016/j.str.2016.10.005."
             <<std::endl<<std::endl;
}

