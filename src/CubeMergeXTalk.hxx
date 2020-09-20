#ifndef CubeMergeXTalk_hxx_seen
#define CubeMergeXTalk_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>
#include <CubeReconCluster.hxx>

#include <memory>
#include <vector>
#include <map>
#include <set>

namespace Cube {
    class MergeXTalk;
};

/// This takes a algorithm result with a TReconObjectContainer with tracks and
/// clusters, and merges any clusters with hits that are consistent with cross
/// talk (or very short delta-rays) into the tracks.  Any tracks that are
/// changed will be refit.  The remaining hits are reclustered so that the
/// clusters do not contain any hits that are also part of a track.
class Cube::MergeXTalk
    : public Cube::Algorithm {
public:

    // A vector of sets of hits associated with each node of a track.
    typedef std::set<Cube::Handle<Cube::Hit>> HitSet;
    // A shared pointer for the set of hits in a node.
    typedef std::shared_ptr<HitSet> NodeHits;
    // All the nodes for one track.
    typedef std::vector< NodeHits > TrackNodeHits;
    // The nodes of all the tracks being analyzed.  One entry per track.
    typedef std::vector< TrackNodeHits > AllTracks;
    // A map from a hit to all of the Nodes that contain the hit.  A hit can
    // be in two or more nodes when tracks overlap (e.g. at a vertex).
    typedef std::map<Cube::Handle<Cube::Hit>, std::set<NodeHits> > AllHits;

    MergeXTalk();
    virtual ~MergeXTalk();

    /// Apply the algorithm.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

private:

    // Count the number of cubes in allHits that are neighboring to hit.
    int CountHitNeighbors(AllHits& allHits, Cube::Handle<Cube::Hit>& hit);

    // Count the number of cubes in a cluster that have a neighbor.
    int CountClusterNeighbors(AllHits& allHits,
                              Cube::Handle<Cube::ReconCluster>& cluster);

    // Find the best neighboring hit to be cross talk.
    Cube::Handle<Cube::Hit> FindBestNeighbor(AllHits& allHits,
                                             Cube::Handle<Cube::Hit>& hit);


};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
