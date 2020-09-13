#ifndef CubeReconShower_hxx_seen
#define CubeReconShower_hxx_seen

#include "CubeHandle.hxx"
#include "CubeReconObject.hxx"
#include "CubeShowerState.hxx"

namespace Cube {
    class ReconShower;
}

/// A representation of an energy deposition starting at a position and
/// falling within a cone.  This is described by the total amount of energy
/// deposited, the starting position of the cone, the time of the deposit, the
/// direction of the cone axis, and the opening angle of the cone (deposit,
/// position, time, direction, opening angle).
///
/// The Cube::ReconShower class is intended to describe the geometry of the
/// energy deposition in a detector, and not make the association with a
/// particular particle identification.
class Cube::ReconShower: public Cube::ReconObject {
public:
    ReconShower();

    /// copy constructor
    ReconShower(const Cube::ReconShower& shower);

    virtual ~ReconShower();

    /// Return a handle to the state.
    Cube::Handle<Cube::ShowerState> GetState() const {
        return GetReconState();
    }

    /// Get the energy deposited in the shower.
    double GetEDeposit() const;

    /// Get the shower starting position.
    TLorentzVector GetPosition() const;

    /// Get the shower starting position uncertainty.
    TLorentzVector GetPositionVariance() const;

    /// Get the number of (non-free) spacial dimensions
    int GetDimensions() const;

    /// Check if this shower has X information.
    bool IsXShower() const;

    /// Check if this shower has Y information.
    bool IsYShower() const;

    /// Check if this shower has Z information.
    bool IsZShower() const;

    /// Get the shower direction.
    TVector3 GetDirection() const;

    /// Get the shower opening angle
    double GetConeAngle() const;

    ClassDef(ReconShower,1);
};
#endif
