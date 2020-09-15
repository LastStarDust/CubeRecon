#ifndef CubePCATrackFit_hxx_seen
#define CubePCATrackFit_hxx_seen

#include "CubeTrackFitBase.hxx"

#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class PCATrackFit;
};

/// Fit a track using a principal component analysis of the track node
/// objects.  This has the effect of robustly fitting a straight line to the
/// track, but is done in a way that does not have a preferred fit direction.
/// The node states are taken by projecting them on the the best fit line.
/// This takes an input track where all of the nodes are filled with objects,
/// and the nodes are in the correct order.  See the Cube::TrackFitBase
/// class for more detailed API documentation.
///
/// Note: This class is a good starting place for more complex fitters since
/// most of what it is doing is to correctly fill the track information.
class Cube::PCATrackFit : public Cube::TrackFitBase {
public:
    PCATrackFit();
    virtual ~PCATrackFit();

    /// Fit the skeleton of a track using a PCA of the node object positions.
    virtual Cube::Handle<Cube::ReconTrack>
    Apply(Cube::Handle<Cube::ReconTrack>& input);
};

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
