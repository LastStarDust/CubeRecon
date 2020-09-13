#include "CubeShowerState.hxx"

ClassImp(Cube::ShowerState);

Cube::ShowerState::ShowerState() {
    ENERGY_DEPOSIT_STATE_DEFINITION;
    POSITION_STATE_DEFINITION;
    DIRECTION_STATE_DEFINITION;
    CONE_STATE_DEFINITION;
    Init();
}

Cube::ShowerState::~ShowerState() {}

Cube::ShowerState::ShowerState(const Cube::ShowerState& init)
    : Cube::ReconState(init) {
    ENERGY_DEPOSIT_STATE_DEFINITION;
    POSITION_STATE_DEFINITION;
    DIRECTION_STATE_DEFINITION;
    CONE_STATE_DEFINITION;
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

Cube::ShowerState&
Cube::ShowerState::operator=(const Cube::ShowerState& rhs) {
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
