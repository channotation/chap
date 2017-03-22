#ifndef ANALYSYS_DATA_LONG_FORMAT_PLOT_MODULE_HPP
#define ANALYSYS_DATA_LONG_FORMAT_PLOT_MODULE_HPP

#include <fstream> 

#include <gromacs/analysisdata/modules/plot.h>


/*
 *
 */
class AnalysisDataLongFormatPlotModule : public gmx::AnalysisDataModuleSerial
{
    public:

        // constructor:
        AnalysisDataLongFormatPlotModule();                                               
        explicit AnalysisDataLongFormatPlotModule(int i); 

        // virtual methods from base class:
        virtual int flags() const;
        virtual void pointsAdded(const gmx::AnalysisDataPointSetRef &points);
        virtual void dataStarted(gmx::AbstractAnalysisData *data);
        virtual void frameStarted(const gmx::AnalysisDataFrameHeader &frame);
        virtual void frameFinished(const gmx::AnalysisDataFrameHeader &frame);
        virtual void dataFinished();

        // setter methods for parameters:
        void setFileName(char *name){fileName_ = name;};
        void setPrecision(int precision){precision_ = precision;};


    private:
    
        // filestream object for output:
        std::fstream file_;

        // internal parameters:
        char *fileName_;
        int precision_;
};

typedef std::shared_ptr<AnalysisDataLongFormatPlotModule> AnalysisDataLongFormatPlotModulePointer;



#endif

