#include <iostream>
#include <fstream>
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
AnalysisDataLongFormatPlotModule::AnalysisDataLongFormatPlotModule(int i)
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
AnalysisDataLongFormatPlotModule::pointsAdded(const gmx::AnalysisDataPointSetRef &points)      
{
    // check that file is open:
    if( file_.is_open() == true ) 
    {
        // output formatiing:
        file_.precision(precision_);

        // write time stamp:
        file_<<points.header().x()<<"\t";
        
        // loop over columns in data set:
        for(unsigned int i = 0; i < points.columnCount() - 1; i++)
        {
            file_<<points.values()[i].value()<<"\t";
        }
        file_<<points.values()[points.lastColumn()].value();

        // line break;
        file_<<std::endl;
    }
} 


/*
 *
 */
void 
AnalysisDataLongFormatPlotModule::dataStarted(gmx::AbstractAnalysisData *data)
{
    std::cout<<"Data started!"<<std::endl;

    // open file stream:
    file_.open(fileName_, std::fstream::out);
}


/*
 *
 */
void 
AnalysisDataLongFormatPlotModule::frameStarted(const gmx::AnalysisDataFrameHeader &frame)
{
    std::cout<<"  Frame started."<<std::endl;
}


/*
 *
 */
void 
AnalysisDataLongFormatPlotModule::frameFinished(const gmx::AnalysisDataFrameHeader &frame)
{
    std::cout<<"  Frame finished."<<std::endl;
}


/*
 *
 */
void 
AnalysisDataLongFormatPlotModule::dataFinished()
{
    std::cout<<"Data finished!"<<std::endl;

    // close file stream:
    file_.close();
}



