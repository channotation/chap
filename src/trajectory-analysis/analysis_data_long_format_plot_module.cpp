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
 * Opens filestream for output writing.
 */
void 
AnalysisDataLongFormatPlotModule::dataStarted(gmx::AbstractAnalysisData *data)
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
AnalysisDataLongFormatPlotModule::frameStarted(const gmx::AnalysisDataFrameHeader &frame)
{

}


/*
 * Nothing to do here.
 */
void 
AnalysisDataLongFormatPlotModule::frameFinished(const gmx::AnalysisDataFrameHeader &frame)
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



