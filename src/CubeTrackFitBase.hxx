#ifndef CubeTrackFitBase_hxx_seen
#define CubeTrackFitBase_hxx_seen

#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class TrackFitBase;
};

/// A base class for track fitting. Classes derived from this class implement
/// a specific track fitting algorithm.
///
/// This fits a skeleton of a track.  The track is expected to have nodes
/// constructed with an Cube::TrackState object and an object derived from
/// Cube::ReconObject.  The nodes must be in order from one end of the track to
/// the other.  The input track is expected to be modified by the fitter so
/// that the result handle will be equal to the input handle.  However, this
/// is not guarranteed.  The result track may be a different object than the
/// input track.  If the fit fails, this returns a NULL handle.
///
/// Most code should be using Cube::TrackFit which will choose the best fitter
/// to use in each circumstance.  How to use these fitting classes:
///
/// \code
/// Cube::TrackFit trackFit;
/// Cube::Handle<Cube::ReconTrack> fittedTrack = trackFit(inputTrack);
/// if (fittedTrack) std::cout << "fit was successful" << std::endl;
/// if (!fittedTrack) std::cout << "fit failed" << std::endl;
/// \endcode
///
/// If the fit fails then the returned handle will be empty.  The track fit
/// can also be applied using the Apply method which will be more convenient
/// when the fitter is referenced by a pointer.
///
/// \code
/// std::unique_ptr<Cube::TrackFitBase> trackFit(new Cube::TrackFit);
/// Cube::Handle<TReconTrack> fittedTrack = trackFit->Apply(inputTrack);
/// \endcode
///
/// \warning The input track is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guarranteed.  The result track may be a different object than the input
/// track.
class Cube::TrackFitBase {
public:
    virtual ~TrackFitBase() {}

    /// Fit the skeleton of a track.  This returns a NULL handle if the fit
    /// fails.  The input track is expected to be modified during the fit.
    /// See Cube::TrackFitBase::Apply() for details.
    Cube::Handle<Cube::ReconTrack>
    operator ()(Cube::Handle<Cube::ReconTrack>& input) {return Apply(input);}

    /// Fit the skeleton of a track.  The track is expected to have nodes
    /// constructed with an Cube::TrackState object and an object derived from
    /// Cube::ReconObject (usually a TReconCluster).  The nodes must be in order
    /// from one end of the track to the other.  The input track is expected
    /// to be modified by the fitter so that the result handle will be equal
    /// to the input handle.  However, this is not guarranteed.  The result
    /// track may be a different object than the input track.  If the fit
    /// fails, this returns a NULL handle.
    virtual Cube::Handle<Cube::ReconTrack>
    Apply(Cube::Handle<Cube::ReconTrack>& input) = 0;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
