#ifndef CHAP_TOPOLOGY_PROVIDER
#define CHAP_TOPOLOGY_PROVIDER

#include <gromacs/selection/selectionoptionbehavior.h>
#include "gromacs/topology/topology.h" 

/*
 *
 */
class ChapTopologyProvider : public gmx::ITopologyProvider
{
    public:

        virtual int getAtomCount()
        {
            return 0;
        }


        virtual t_topology *getTopology(bool required)
        {
            
        }
};

#endif

