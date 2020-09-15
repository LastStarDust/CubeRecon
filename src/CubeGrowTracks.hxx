#ifndef CubeGrowTracks_hxx_seen
#define CubeGrowTracks_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>
#include <CubeReconTrack.hxx>
#include <CubeTrackState.hxx>
#include <CubeReconCluster.hxx>

namespace Cube {
    class GrowTracks;
};

/// This takes a algorithm result with a TReconObjectContainer of with track
/// segments and merge them into longer tracks.  The track segments should be
/// fit prior to running this algorithm.  Only the states for the first and
/// last nodes are used.
class Cube::GrowTracks
    : public Cube::Algorithm {
public:
    GrowTracks();
    virtual ~GrowTracks();

    /// Apply the algorithm.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// An enum specifiying how the tracks are oriented.
    typedef enum {kFrontFront, kFrontBack, kBackFront, kBackBack, kNotClose}
        Orientation;

private:
    /// Return how the tracks are oriented relative to each other.
    Orientation TrackOrientation(Cube::Handle<Cube::ReconTrack> t1,
                                 Cube::Handle<Cube::ReconTrack> t2);

    // Return the states that are matchd for tracks t1 and t2.  The directions
    // need to be reversed if the reverseLeading and reverseFollowing are
    // negative.
    Orientation MatchedStates(Cube::Handle<Cube::ReconTrack> t1,
                              Cube::Handle<Cube::ReconTrack> t2,
                              Cube::Handle<Cube::TrackState>& leading,
                              double& reverseLeading,
                              Cube::Handle<Cube::TrackState>& following,
                              double& reverseFollowing);

    /// Return the distance to the point of closest approach.  A positive
    /// value means the position is in the direction that dir is pointing.
    double DistToApproach(const TVector3& dest,
                          const TVector3& src,
                          const TVector3& dir);

    /// Return the distance of closest approach.  This will always be positive.
    double ApproachDist(const TVector3& dest,
                        const TVector3& src,
                        const TVector3& dir);

    /// Return a goodness for how well the tracks are matched.  The goodness
    /// is based on the likelihood that the ends of the tracks line up.  A
    /// higher value is a worse match.
    double MatchGoodness(Cube::Handle<Cube::ReconTrack> t1,
                         Cube::Handle<Cube::ReconTrack> t2);

    /// Take two tracks and merge the clusters in the right order to build a
    /// third track which is returned.  If the tracks can't be merged this
    /// will return an empty handle.
    Cube::Handle<Cube::ReconTrack>
    MergeTracks(Cube::Handle<Cube::ReconTrack> t1,
                Cube::Handle<Cube::ReconTrack> t2);

    /// Return the chi2 for the three clusters falling in a line.  The chi2
    /// will have 3 d.o.f. since there are 9 measurements, and 6 parameters.
    double ThreeInLine(Cube::Handle<Cube::ReconCluster> a,
                       Cube::Handle<Cube::ReconCluster> b,
                       Cube::Handle<Cube::ReconCluster> c);

    /// The cut value for the maximum distance between clusters at the end of
    /// a track to be in the same track.
    double fMergeDistanceCut;

    /// The cut value for the minimum value of the distance from the tail of
    /// the front track to the head of the following track.  The apriori value
    /// is about negative one centimeter.
    double fMinimumFollowDistanceCut;

    /// Tracks that have a kink with a cosine of less than this are not even
    /// considered for merging.
    double fFollowCosineCut;

    /// The cut value for the chi2 for clusters to be considered as from
    /// the same line.
    double fGoodnessCut;

    /// Tracks have directions more closely matched than this are not
    /// rejected.  They will also need to be close together.
    double fMatchedCosineCut;

    /// Tracks that have a closets approach of less than this are not
    /// rejected.  They will also need to have matched directions.
    double fMatchedPositionCut;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
