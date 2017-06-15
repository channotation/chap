#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/exceptions.h"

#include "commandline/chap_topology_information.hpp"



/*
 *
 */
ChapTopologyInformation::ChapTopologyInformation()
    : top_(NULL)
    , bTop_(false)
    , xtop_(NULL)
    , ePBC_(-1)
{
    clear_mat(boxtop_);
}


/*
 *
 */
ChapTopologyInformation::~ChapTopologyInformation()
{
    if( top_ )
    {
        done_top(top_);
        sfree(top_);
    }
    sfree(xtop_);
}


/*
 *
 */
void
ChapTopologyInformation::getTopologyConf(rvec **x, matrix box) const
{
    if( box )
    {
        copy_mat(const_cast<rvec *>(boxtop_), box);
    }
    if( x )
    {
        if( !xtop_ )
        {
            *x = NULL;
            GMX_THROW(gmx::APIError("Topology coordinates requested without setting efUseTopX"));
        }
        *x = xtop_;
    }
}

