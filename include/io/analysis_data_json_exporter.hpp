#ifndef ANALYSIS_DATA_JSON_EXPORTER
#define ANSLYSIS_DATA_JSON_EXPORTER

#include <string>
#include <vector>

#include <gromacs/analysisdata/datamodule.h>

#include "rapidjson/document.h"


/*
 *
 */
class AnalysisDataJsonExporter : public gmx::AnalysisDataModuleParallel 
{
    public:

        // constructor and destructor:
        AnalysisDataJsonExporter(){};
        ~AnalysisDataJsonExporter(){};

        // implementation of interface as specified by base class:
        virtual int flags() const;
        virtual bool parallelDataStarted(gmx::AbstractAnalysisData *data,
                                         const gmx::AnalysisDataParallelOptions &options);
        virtual void frameStarted(const gmx::AnalysisDataFrameHeader &frame);
        virtual void pointsAdded(const gmx::AnalysisDataPointSetRef &points);
        virtual void frameFinished(const gmx::AnalysisDataFrameHeader &frame);
        virtual void frameFinishedSerial(int index);
        virtual void dataFinished();

        // setter methods for parameters:
        void setDataSetNames(std::vector<std::string> dataSetNames);

    private:

        // json document to carry data:
        rapidjson::Document json_;

        // internal parameters:
        std::vector<std::string> dataSetNames_;
};


/*
 *
 */
typedef std::shared_ptr<AnalysisDataJsonExporter> AnalysisDataJsonExporterPointer;

#endif

