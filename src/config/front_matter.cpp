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

