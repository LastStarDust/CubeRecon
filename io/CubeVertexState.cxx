#include "CubeVertexState.hxx"

///////////////////////////////////////////////////////
ClassImp(Cube::VertexState);

Cube::VertexState::VertexState() {

    POSITION_STATE_DEFINITION;

    Init();
}

Cube::VertexState::~VertexState() {}

Cube::VertexState::VertexState(const Cube::VertexState& init)
    : Cube::ReconState(init) {

    POSITION_STATE_DEFINITION;

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

Cube::VertexState&
Cube::VertexState::operator=(const Cube::VertexState& rhs) {
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
