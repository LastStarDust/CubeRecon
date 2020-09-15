#include <CubeCreateTrack.hxx>
#include <CubeReconTrack.hxx>
#include <CubeTrackState.hxx>
#include <CubeReconNode.hxx>

Cube::Handle<Cube::ReconTrack>
Cube::CombineTracks(const Cube::ReconTrack& t1,
                    const Cube::ReconTrack& t2) {

    Cube::Handle<Cube::TrackState> s1Front = t1.GetNodes().front()->GetState();
    Cube::Handle<Cube::TrackState> s1Back = t1.GetNodes().back()->GetState();
    Cube::Handle<Cube::TrackState> s2Front = t2.GetNodes().front()->GetState();
    Cube::Handle<Cube::TrackState> s2Back = t2.GetNodes().back()->GetState();


    double lenS1 = (s1Front->GetPosition().Vect()
                    - s1Back->GetPosition().Vect()).Mag();
    double lenS2 = (s1Front->GetPosition().Vect()
                    - s1Back->GetPosition().Vect()).Mag();
    double distS1FrontToS2Front = (s1Front->GetPosition().Vect()
                               - s2Front->GetPosition().Vect()).Mag();
    double distS1FrontToS2Back = (s1Front->GetPosition().Vect()
                              - s2Back->GetPosition().Vect()).Mag();
    double distS1BackToS2Front = (s1Back->GetPosition().Vect()
                              - s2Front->GetPosition().Vect()).Mag();
    double distS1BackToS2Back = (s1Back->GetPosition().Vect()
                             - s2Back->GetPosition().Vect()).Mag();

    typedef enum {kFrontToFront,
                  kFrontToBack,
                  kBackToFront,
                  kBackToBack} Match;
    Match match;

    if (distS1BackToS2Front < distS1FrontToS2Front
             && distS1BackToS2Front < distS1FrontToS2Back
             && distS1BackToS2Front < distS1BackToS2Back) {
        match = kBackToFront;
    }
    else if (distS1BackToS2Back < distS1FrontToS2Front
             && distS1BackToS2Back < distS1FrontToS2Back
             && distS1BackToS2Back < distS1BackToS2Front) {
        match = kBackToBack;
    }
    else if (distS1FrontToS2Front < distS1FrontToS2Back
        && distS1FrontToS2Front < distS1BackToS2Front
        && distS1FrontToS2Front < distS1BackToS2Back) {
        match = kFrontToFront;
    }
    else if (distS1FrontToS2Back < distS1FrontToS2Front
             && distS1FrontToS2Back < distS1BackToS2Front
             && distS1FrontToS2Back < distS1BackToS2Back) {
        match = kFrontToBack;
    }
    else {
        CUBE_ERROR << "Inconceivable!" << std::endl;
    }

    std::vector< Cube::Handle<Cube::ReconCluster> > clusters;

    // Add the clusters from the first track.
    if (match == kBackToBack || match == kBackToFront) {
        // The front track (track1) is already going in the correct direction.
        for (Cube::ReconNodeContainer::const_iterator
                 n = t1.GetNodes().begin();
             n != t1.GetNodes().end(); ++n) {
            clusters.push_back((*n)->GetObject());
        }
    }
    else {
        // The front track (track1) is going in the wrong direction, so it
        // needs to be reversed.
        for (Cube::ReconNodeContainer::const_reverse_iterator
                 n = t1.GetNodes().rbegin();
             n != t1.GetNodes().rend(); ++n) {
            clusters.push_back((*n)->GetObject());
        }
    }

    // Add the clusters from the second track.
    if (match == kFrontToFront || match == kBackToFront) {
        // The back track (track2) is going in the correct direction.
        for (Cube::ReconNodeContainer::const_iterator
                 n = t2.GetNodes().begin();
             n != t2.GetNodes().end(); ++n) {
            clusters.push_back((*n)->GetObject());
        }
    }
    else {
        // The back track (track2) is going in the wrong direction, so it
        // needs to be reversed.
        for (Cube::ReconNodeContainer::const_reverse_iterator
                 n = t2.GetNodes().rbegin();
             n != t2.GetNodes().rend(); ++n) {
            clusters.push_back((*n)->GetObject());
        }
    }

    // Construct the new track.
    Cube::Handle<Cube::ReconTrack> result;
    result = Cube::CreateTrackFromClusters(
        "combined", clusters.begin(), clusters.end());

    // Copy the very front and back states from the old tracks so that the new
    // track has the best estimate of the states at the ends.  This is not
    // equivalent to refitting the track, but may help with efficiency since
    // the front and back nodes will have the best current estimated state.

    if (match == kBackToFront || match == kBackToBack) {
        // Set the front state from the front of track1
        result->GetNodes().front()->GetState()
            = t1.GetNodes().front()->GetState();
    }
    else {
        // Set the front state from the back of track1 and reverse the
        // direction.
        result->GetNodes().front()->GetState()
            = t1.GetNodes().back()->GetState();
        Cube::Handle<Cube::TrackState> state
            = result->GetNodes().front()->GetState();
        state->SetDirection(-state->GetDirection());
    }

    if (match == kBackToFront || match == kBackToBack) {
        // Set the back state from the back of track2
        result->GetNodes().back()->GetState()
            = t2.GetNodes().back()->GetState();
    }
    else {
        // Set the back state from the front of track2 and reverse the
        // direction.
        result->GetNodes().back()->GetState()
            = t2.GetNodes().front()->GetState();
        Cube::Handle<Cube::TrackState> state
            = result->GetNodes().back()->GetState();
        state->SetDirection(-state->GetDirection());
    }

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
