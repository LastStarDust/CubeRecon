#ifndef CubeShowerState_hxx_seen
#define CubeShowerState_hxx_seen

#include "CubeReconState.hxx"

namespace Cube {
    class ShowerState;
}

/// A state holding the parameters associated with a TReconShower.
class Cube::ShowerState:
    public Cube::ReconState {
public:
    ShowerState();
    ShowerState(const ShowerState& init);
    virtual ~ShowerState();
    virtual ShowerState& operator=(const ShowerState& rhs);

    ENERGY_DEPOSIT_STATE_DECLARATION;
    POSITION_STATE_DECLARATION;
    DIRECTION_STATE_DECLARATION;
    CONE_STATE_DECLARATION;

    ENERGY_DEPOSIT_STATE_PRIVATE;
    POSITION_STATE_PRIVATE;
    DIRECTION_STATE_PRIVATE;
    CONE_STATE_PRIVATE;

    ClassDef(ShowerState,1);
};
#endif
