#ifndef CubeVertexState_hxx_seen
#define CubeVertexState_hxx_seen

#include "CubeReconState.hxx"

namespace Cube {
    class VertexState;
}

/// A state holding parameters associated with a TReconVertex.
class Cube::VertexState:
    public Cube::ReconState {
public:
    VertexState();
    virtual ~VertexState();
    VertexState(const VertexState& init);
    virtual VertexState& operator=(const VertexState& rhs);

    POSITION_STATE_DECLARATION;

    POSITION_STATE_PRIVATE;

    ClassDef(VertexState,1);
};
#endif
