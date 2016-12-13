#ifndef PATH_FINDING_MODULE_HPP
#define PATH_FINDING_MODULE_HPP

#include <vector>

#include <gromacs/trajectoryanalysis.h>

class PathFindingModule
{
	public:
		
		PathFindingModule(gmx::RVec initProbePos,
                          gmx::RVec chanDirVec,
                          gmx::AnalysisNeighborhoodSearch nbSearch,
                          std::vector<real> vdwRadii);
		~PathFindingModule();

        void findPath();

        int getCulprit(){return culprit;};

        std::vector<gmx::RVec> getPath(){return(path_);};
        std::vector<real> getRadii(){return(radii_);};

//	private:

        int culprit;

        real stepLength_;

        std::vector<real> vdwRadii_;

        gmx::RVec initProbePos_;
        gmx::RVec prevProbePos_;
        gmx::RVec crntProbePos_;

        gmx::RVec chanDirVec_;                  // channel direction vector
        gmx::RVec orthVecU_;                    // first orthogonal vector
        gmx::RVec orthVecW_;                    // second orthogonal vector

        const gmx::AnalysisNeighborhoodSearch nbSearch_;

        std::vector<gmx::RVec> path_;
        std::vector<real> radii_;

        void optimiseInitial();
        void marchAndOptimise(gmx::RVec initPos, bool forward);

        real findMinimumDistance(real *state);

        gmx::RVec optimToConfig(real *state);


};





#endif

