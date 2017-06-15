#ifndef CHAP_TOPOLOGY_INFORMATION_HPP
#define CHAP_TOPOLOGY_INFORMATION_HPP

#include <gromacs/math/vectypes.h>
#include <gromacs/topology/topology.h>


/*
 * 
 */
class ChapTopologyInformation
{
    public:

        ChapTopologyInformation();
        ~ChapTopologyInformation();

        bool hasTopology() const { return top_ != NULL; };
        bool hasFullTopology() { return bTop_; };

        t_topology *topology() const { return top_; };
        int ePBC() const { return ePBC_; };

        void getTopologyConf(rvec **x, matrix box) const;



    private:


        t_topology *top_;
        bool bTop_;
        rvec *xtop_;
        matrix boxtop_;
        int ePBC_;
        
        friend class ChapTrajAnaRunnerCommon;
};


#endif

