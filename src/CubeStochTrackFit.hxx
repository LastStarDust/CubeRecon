#ifndef CubeStochTrackFit_hxx_seen
#define CubeStochTrackFit_hxx_seen

#include "CubeTrackFitBase.hxx"

#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class StochTrackFit;
};

/// A stochastic track fit done with a Sequential Importance Resampling (SIR)
/// Particle Filter using a Gaussian uncertainty around each hit in the track.
/// The SIR filter is implemented using the SimpleSIR template.  This takes an
/// input track where all of the nodes are filled with objects, and the nodes
/// are in the correct order.  See the Cube::TrackFitBase class for more
/// detailed API documentation.
///
/// Note: This class is a good starting place for more complex SIR filters
/// since most of what it is doing is to correctly fill the track information.
class Cube::StochTrackFit : public Cube::TrackFitBase {
public:

    /// Construct the track fitter (usually done inside TrackFit).  This
    /// takes an optional integer which will determine the number of samples
    /// to be used to approximate the PDF.
    explicit StochTrackFit(int nSamples = 1000);
    virtual ~StochTrackFit();

    /// Fit the skeleton of a track.
    virtual Cube::Handle<Cube::ReconTrack>
    Apply(Cube::Handle<Cube::ReconTrack>& input);

    /// Set the number of samples used to approximate the PDF.  The number of
    /// samples will be used in the next call to Apply.
    void SetSampleCount(int s) {fSampleCount = s;}

    /// Get the number of samples that will be used in the next call to Apply.
    int GetSampleCount() const {return fSampleCount;}

private:
    // The number of samples in the sample vector that is used to describe the
    // PDF.
    int fSampleCount;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
