#ifndef ANALYSYS_DATA_PDB_PLOT_MODULE_HPP
#define ANALYSYS_DATA_PDB_PLOT_MODULE_HPP

#include <fstream>
#include <vector>

#include <gromacs/analysisdata/modules/plot.h>


/*
 *
 */
class AnalysisDataPdbPlotModule : public gmx::AnalysisDataModuleSerial
{
    public:

        // constructor:
        AnalysisDataPdbPlotModule();                                               
        explicit AnalysisDataPdbPlotModule(int i); 

        // virtual methods from base class:
        virtual int flags() const;
        virtual void pointsAdded(const gmx::AnalysisDataPointSetRef &points);
        virtual void dataStarted(gmx::AbstractAnalysisData *data);
        virtual void frameStarted(const gmx::AnalysisDataFrameHeader &frame);
        virtual void frameFinished(const gmx::AnalysisDataFrameHeader &frame);
        virtual void dataFinished();

        // setter methods for parameters:
        void setFileName(const char *name){fileName_ = name;};
        void setPrecision(int precision){precision_ = precision;};
        void setHeader(std::vector<char*> header){header_ = header;};


    private:
    
        // filestream object for output:
        std::fstream file_;

        // internal parameters:
        const char *fileName_;
        int precision_;
        std::vector<char*> header_;
};

typedef std::shared_ptr<AnalysisDataPdbPlotModule> AnalysisDataPdbPlotModulePointer;



#endif

