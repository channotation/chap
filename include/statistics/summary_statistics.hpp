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


#ifndef SUMMARY_STATISTICS_HPP
#define SUMMARY_STATISTICS_HPP

#include <vector>

#include "gromacs/utility/real.h"


/*!
 * \brief Collects summary statistics of a scalar variable without having to
 * hold an entire dataset in memory.
 *
 * This class computes the minimum, maximum, mean, variance, standard 
 * deviation, and number of samples of a scalar-valued dataset in an incremental
 * fashion. It exposes an update() method which is used to update the summary
 * statistics with new values, but does not store the values themselves.
 *
 * For the minimum and maximum this is accomplished by a trivial comparison to
 * the currently stored minimum and maximum. For mean, variance, and standard
 * deviation a numerically stable algorithm due to Welford (1962) is used. 
 *
 * Note that while standard deviation and variance are strictly speaking 
 * undefined for less then two data points, this class will return a value of
 * zero in this case to simplify data handling in the context of JSON.
 *
 * Note furthermore that this class will always skip infinite values and in 
 * this case no variable is updated (the counter also stays at zero). As a
 * precaution for situation where no (finite) values have been passed to 
 * the update() method, all returned statistics are capped at the minimum and 
 * maximum real number, i.e. no infinity is ever returned to allow 
 * compatibility with JSON.
 */
class SummaryStatistics
{
    public:

        // constructor and destructor:
        SummaryStatistics();        

        // getter methods:
        real min() const;
        real max() const;
        real mean() const;
        real var() const;
        real sd() const;
        int num() const;

        // updating method:
        void update(
                const real newValue);
        static void updateMultiple(
                std::vector<SummaryStatistics> &stat,
                const std::vector<real> &newValues);

        // manipulation methods:
        void shift(
                const real shift);

    private:

        // summary statistics updated by this class:
        real min_;
        real max_;
        real mean_;
        real sumSquaredMeanDiff_;
        int num_;

        // internal auxiliary functions:
        inline real varFromSumSquaredMeanDiff() const;
        inline real mendInfinity(real value) const;
};


#endif

