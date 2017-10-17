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

