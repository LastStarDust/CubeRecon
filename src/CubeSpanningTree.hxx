#ifndef CubeSpanningTree_hxx_seen
#define CubeSpanningTree_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class SpanningTree;
};

/// Accept a input result with clusters and break the clusters up into smaller
/// ones based on the junctions of a minimal spanning tree.  This produces an
/// object container with clusters.  The clusters have the hits ordered
/// according to the order in the minimal spanning tree.  The first and last
/// hits each cluster correspond to the hits at the branching points in the
/// tree.  The hit at a branching point is in each cluster associated with the
/// branching point.  In "graph" language, the clusters are (kind-of-like) the
/// edges, and the end points of the clusters (the duplicated hits) are
/// (kind-of-like) the vertices.
class Cube::SpanningTree: public Cube::Algorithm {
public:
    SpanningTree();
    virtual ~SpanningTree() {}

    /// The routine that does the actual work.
    Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    void SetDistanceType(int i) {fDistanceType = i;}

    void SetOversizeCut(int i) {fOversizeCut = i;}

private:

    // Control how the edge weight is calculated.
    int fDistanceType;

    // The spanning tree implementation can be pretty slow, so protect against
    // really large number of hits [the current simple implementation is
    // O(N^3)!!!].  This keeps the code from hanging on large events.  It
    // should be as large as possible, while preventing the reconstruction
    // from taking several 10's of minutes for the big clusters.
    int fOversizeCut;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
