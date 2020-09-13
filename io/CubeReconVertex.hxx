#ifndef CubeReconVertex_hxx_seen
#define CubeReconVertex_hxx_seen

#include "CubeHandle.hxx"
#include "CubeReconObject.hxx"
#include "CubeVertexState.hxx"

namespace Cube {
    class ReconVertex;
}

/// Define a vertex location within the detector.
class Cube::ReconVertex: public Cube::ReconObject {
public:
    ReconVertex();
    virtual ~ReconVertex();

    /// Return a handle to the state.
    Cube::Handle<Cube::VertexState> GetState() const {
        return GetReconState();
    }

    /// Get the vertex position.
    TLorentzVector GetPosition() const;

    /// Get the track starting position uncertainty.
    TLorentzVector GetPositionVariance() const;

    /// Get the number of (non-free) spacial dimensions
    int GetDimensions() const;

    /// Check if this vertex has X information.
    bool IsXVertex() const;

    /// Check if this vertex has Y information.
    bool IsYVertex() const;

    /// Check if this vertex has Z information.
    bool IsZVertex() const;

    ClassDef(ReconVertex,1);
};
#endif
