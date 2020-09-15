#ifndef CubeTrackFit_hxx_seen
#define CubeTrackFit_hxx_seen

#include "CubeTrackFitBase.hxx"

#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class TrackFit;
    class PCATrackFit;
    class StochTrackFit;
};

/// A class to fit the skeleton of a track.  The track is expected to have
/// nodes constructed with a Cube::TTrackState object and an object derived
/// from Cube::ReconObject.  The nodes must be in order from one end of the
/// track to the other.  The input track is expected to be modified by the
/// fitter so that the result handle will be equal to the input handle.
/// However, this is not guarranteed.  The result track may be a different
/// object than the input track.  If the fit fails, this returns a NULL
/// handle.
///
/// This is a wrapper around other track fitting classes (all derived from
/// Cube::TrackFitBase that chooses the correct fitter to be applied.  The
/// Cube::TrackFit class is the "main" class serving as a switch yard to
/// determine the best fitter for each type of track.  Notice that this is
/// fitting ReconTrack objects, not ReconPID objects.  ReconPID objects
/// must be fit with a different class.
///
/// Most code should be using Cube::TrackFit which will choose the best fitter
/// to use in each circumstance.  How to use these fitting classes:
///
/// \code
/// Cube::TrackFit trackFit;
/// Cube::Handle<ReconTrack> fittedTrack = trackFit(inputTrack);
/// if (fittedTrack) std::cout << "fit was successful" << std::endl;
/// if (!fittedTrack) std::cout << "fit failed" << std::endl;
/// \endcode
///
/// If the fit fails then the returned handle will be empty.  The track fit
/// can also be applied using the Apply method which will be more convenient
/// when the fitter is referenced by a pointer.
///
/// \code
/// std::unique_ptr<Cube::TrackFit> trackFit(new Cube::TrackFit);
/// Cube::Handle<ReconTrack> fittedTrack = trackFit->Apply(inputTrack);
/// \endcode
///
/// \warning The input track is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guarranteed.  The result track may be a different object than the input
/// track.
class Cube::TrackFit : public Cube::TrackFitBase {
public:
    explicit TrackFit();
    virtual ~TrackFit();

    /// Fit the skeleton of a track.  The track is expected to have nodes
    /// constructed with a Cube::TTrackState and an object derived from
    /// Cube::ReconObject.  The nodes must be in order from one end of the
    /// track to the other.  The input track is expected to be modified by the
    /// fitter so that the result handle will be equal to the input handle.
    /// However, this is not guarranteed.  The result track may be a different
    /// object than the input track.  If the fit fails, this returns a NULL
    /// handle.
    virtual Cube::Handle<Cube::ReconTrack>
    Apply(Cube::Handle<Cube::ReconTrack>& input);

    // Return a pointer to the PCA track fitter (it may be NULL).
    Cube::PCATrackFit* GetPCA() {return fPCA;}

    // Return a pointer to the stochastic track fitter (it may be NULL).
    Cube::StochTrackFit* GetStochastic() {return fStochastic;}

private:

    /// A pointer to the stochastic stocastic fitter.  This is only
    /// instantiated if the fitter is used.
    Cube::StochTrackFit* fStochastic;

    /// A pointer to the PCA fitter.  This is only instantiated if the
    /// fitter is used.
    Cube::PCATrackFit* fPCA;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
