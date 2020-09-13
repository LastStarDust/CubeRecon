#ifndef CubeClusterState_hxx_seen
#define CubeClusterState_hxx_seen

#include "CubeReconState.hxx"

namespace Cube {
    class ClusterState;
}

/// A state holding the parameters associated with a TReconCluster.
class Cube::ClusterState:
    public Cube::ReconState
{
public:
    ClusterState();
    virtual ~ClusterState();
    ClusterState(const ClusterState& init);
    virtual ClusterState& operator=(const ClusterState& rhs);

    ENERGY_DEPOSIT_STATE_DECLARATION;
    POSITION_STATE_DECLARATION;

    ENERGY_DEPOSIT_STATE_PRIVATE;
    POSITION_STATE_PRIVATE;

    ClassDef(ClusterState,1);
};
#endif
