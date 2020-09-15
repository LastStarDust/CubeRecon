#include "CubeTrackFit.hxx"
#include "CubePCATrackFit.hxx"
#include "CubeStochTrackFit.hxx"

Cube::TrackFit::TrackFit()
    : fStochastic(NULL), fPCA(NULL) {}
Cube::TrackFit::~TrackFit() {
    if (fStochastic) delete fStochastic;
    if (fPCA) delete fPCA;
}

Cube::Handle<Cube::ReconTrack>
Cube::TrackFit::Apply(Cube::Handle<Cube::ReconTrack>& input) {
    Cube::Handle<Cube::ReconTrack> result;

    // Try to apply the stochastic stocastic fitter.
    if (!fStochastic) {
        fStochastic = new Cube::StochTrackFit;
    }

    result = fStochastic->Apply(input);

    // If theres a successful result, return it. (not very useful without
    // multiple fitters!).
    if (result) return result;

    // Try to apply the PCA fitter.  It should always work, but only fits a
    // straight line.
    if (!fPCA) {
        fPCA = new Cube::PCATrackFit;
    }

    result = fPCA->Apply(input);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
