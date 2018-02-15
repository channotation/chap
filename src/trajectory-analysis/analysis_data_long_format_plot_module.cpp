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


#include <cstring>
#include <iomanip>

#include <gromacs/analysisdata/dataframe.h>

#include "trajectory-analysis/analysis_data_long_format_plot_module.hpp"


/*
 *
 */
AnalysisDataLongFormatPlotModule::AnalysisDataLongFormatPlotModule()                                
    : precision_(5)
{

}                                           


/*
 *
 */
AnalysisDataLongFormatPlotModule::AnalysisDataLongFormatPlotModule(int /*i*/)
    : precision_(5)
{

}                                                                               


/*
 *
 */
int
AnalysisDataLongFormatPlotModule::flags() const
{
    return efAllowMissing | efAllowMulticolumn | efAllowMultipoint
           | efAllowMultipleDataSets; 
}


/*
 *
 */
void                                                                            
AnalysisDataLongFormatPlotModule::pointsAdded(
        const gmx::AnalysisDataPointSetRef &points)      
{
    // check that file is open:
    if( file_.is_open() == true ) 
    {
        // write time stamp:
        file_<<points.header().x()<<"\t";
        
        // loop over columns in data set:
        for(int i = 0; i < points.columnCount() - 1; i++)
        {
            file_<<points.values()[i].value()<<"\t";
        }
        file_<<points.values()[points.lastColumn()].value();

        // line break;
        file_<<std::endl;
    }
} 


/*
 * Opens filestream for output writing.
 */
void 
AnalysisDataLongFormatPlotModule::dataStarted(
        gmx::AbstractAnalysisData* /*data*/)
{
    // check that file name is set:
    if( strlen(fileName_) > 0 )
    {
        // open file stream:
        file_.open(fileName_, std::fstream::out);

        // output formatiing:
        file_.precision(precision_);

        // if set, write header:
        if( header_.size() )
        {
            // TODO: assert that header has as many elements as data columns?
            for(unsigned int i = 0; i < header_.size() - 1; i++)
            {
                file_<<header_[i]<<"\t";
            }
            file_<<header_.back()<<std::endl;
        }
    }
}


/*
 * Nothing to do here.
 */
void 
AnalysisDataLongFormatPlotModule::frameStarted(
        const gmx::AnalysisDataFrameHeader& /*frame*/)
{

}


/*
 * Nothing to do here.
 */
void 
AnalysisDataLongFormatPlotModule::frameFinished(
        const gmx::AnalysisDataFrameHeader& /*frame*/)
{

}


/*
 * Closes the file stream when all data has been written.
 */
void 
AnalysisDataLongFormatPlotModule::dataFinished()
{
    // close file stream:
    file_.close();
}



