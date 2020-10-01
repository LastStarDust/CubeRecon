#include "CubeClusterState.hxx"

///////////////////////////////////////////////////////
ClassImp(Cube::ClusterState);

Cube::ClusterState::ClusterState() {

    POSITION_STATE_DEFINITION;
    ENERGY_DEPOSIT_STATE_DEFINITION;

    Init();
}

Cube::ClusterState::~ClusterState() {}

Cube::ClusterState::ClusterState(const Cube::ClusterState& init) {

    POSITION_STATE_DEFINITION;
    ENERGY_DEPOSIT_STATE_DEFINITION;

    Init();

    for (int i=0; i<GetDimensions(); ++i) {
        SetValue(i,init.GetValue(i));
    }

    for (int i=0; i<GetDimensions(); ++i) {
        for (int j=0; j<GetDimensions(); ++j) {
            SetCovarianceValue(i,j,init.GetCovarianceValue(i,j));
        }
    }

}

Cube::ClusterState& Cube::ClusterState::operator=(
    const Cube::ClusterState& rhs) {
    if (this == &rhs) return *this;

    for (int i=0; i<GetDimensions(); ++i) {
        SetValue(i,rhs.GetValue(i));
    }

    for (int i=0; i<GetDimensions(); ++i) {
        for (int j=0; j<GetDimensions(); ++j) {
            SetCovarianceValue(i,j,rhs.GetCovarianceValue(i,j));
        }
    }

    return *this;
}
