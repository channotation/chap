#ifndef MOCK_PORE_MAKER_HPP
#define MOCK_PORE_MAKER_HPP

#include <vector>

#include <gromacs/math/vec.h>



/*
 *
 */
class MockPoreMaker
{
    public:

        std::vector<gmx::RVec> makePore(real poreLength,
                                        real poreCentreRadius,
                                        real poreVdwRadius,
                                        gmx::RVec poreCentre,
                                        int alongAxis);


    private:

        const real PI_ = std::acos(-1.0);

};

#endif

