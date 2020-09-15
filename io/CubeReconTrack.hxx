#ifndef CubeReconTrack_hxx_seen
#define CubeReconTrack_hxx_seen

#include <TLorentzVector.h>
#include <TVector3.h>

#include "CubeHandle.hxx"
#include "CubeReconObject.hxx"
#include "CubeTrackState.hxx"

namespace Cube {
    class ReconTrack;
}

/// A representation of a curvilinear energy deposit starting a position, and
/// following a path.  This is described by the amount of energy deposited by
/// the entire track, the initial position , the initial time, the initial
/// direction, the estimated curvature, and the initial path width.  There
/// must be a way to represent the path of the energy deposition between the
/// initial and final ends of the deposition.  At each intermediate point, we
/// require a representation of the energy deposit (dEdX), position, time,
/// direction, curvature, and width.  The detector hits are associated
/// with each node along the track.
///
/// The Cube::ReconTrack class is intended to describe the geometry of the
/// energy deposition in a detector, and not make the association with a
/// particular particle identification.
class Cube::ReconTrack: public Cube::ReconObject {
public:
    ReconTrack();

    /// copy constructor
    ReconTrack(const Cube::ReconTrack& track);

    virtual ~ReconTrack();

    /// Return a handle to the state.  The state at the front end of the track
    /// is not necessarily the same as the state at the first node.  For
    /// instance, the first node may have a finite extent, and the starting
    /// state of the track will be the estimated state at the edge of the node
    /// extent.  The state of the node will be the state of the track at the
    /// center of the node.
    Cube::Handle<Cube::TrackState> GetState() const {
        return GetReconState();
    }

    /// Return a handle to the state at the front end of the track.  The front
    /// and back of the track are defined by the order of the nodes.  See
    /// GetState() for more details.
    Cube::Handle<Cube::TrackState> GetFront() const {
        return GetState();
    }

    /// Return a handle to the state at the back end of the track.  The front
    /// and back of the track are defined by the order of the nodes.  The
    /// state at the back end of the track may not be the same as the state of
    /// the last node.  See GetState() for more details.
    Cube::Handle<Cube::TrackState> GetBack() const {
        return Cube::Handle<Cube::TrackState>(fBackState,false);
    }

    /// Get the energy deposited in the track.
    double GetEDeposit() const;

    /// Get the track starting position.
    TLorentzVector GetPosition() const;

    /// Get the track starting position uncertainty.
    TLorentzVector GetPositionVariance() const;

    /// Get the number of (non-free) spacial dimensions
    int GetDimensions() const;

    /// Check if this track has X information.
    bool IsXTrack() const;

    /// Check if this track has Y information.
    bool IsYTrack() const;

    /// Check if this track has Z information.
    bool IsZTrack() const;

    /// Get the track direction.
    TVector3 GetDirection() const;

    /// Get the track curvature.
    TVector3 GetCurvature() const;

    /// Get the width of the track.
    double GetWidth() const;

    /// Reverse the track direction.
    void ReverseTrack();

    /// List the results of particle id.
    virtual void ls(Option_t* opt = "") const;

private:

    /// The state of the track at the back end.
    Cube::TrackState* fBackState;

    ClassDef(ReconTrack,2);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
