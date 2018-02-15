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


#ifndef ANALYSYS_DATA_LONG_FORMAT_PLOT_MODULE_HPP
#define ANALYSYS_DATA_LONG_FORMAT_PLOT_MODULE_HPP

#include <fstream>
#include <string>
#include <vector>

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
        void setFileName(const char *name){fileName_ = name;};
        void setPrecision(int precision){precision_ = precision;};
        void setHeader(std::vector<std::string> header){header_ = header;};


    private:
    
        // filestream object for output:
        std::fstream file_;

        // internal parameters:
        const char *fileName_;
        int precision_;
        std::vector<std::string> header_;
};

typedef std::shared_ptr<AnalysisDataLongFormatPlotModule> AnalysisDataLongFormatPlotModulePointer;



#endif

