#ifndef SUMMARY_STATISTICS_HPP
#define SUMMARY_STATISTICS_HPP

#include <vector>

#include "gromacs/utility/real.h"


/*!
 * \brief Collects summary statistics of a scalar variable without having to
 * hold an entire dataset in memory.
 *
 * This class computes the minimum, maximum, mean, variance, standard 
 * devition, and number of samples of a scalar-valued dataset in an incremental
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
 * precaution for situatuation where no (finite) values have been passed to 
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

