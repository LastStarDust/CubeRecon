#ifndef SFGGrowClusters_hxx_seen
#define SFGGrowClusters_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class GrowClusters;
};

/// Grow smaller clusters into larger track-like clusters.  This take a group
/// of clusters that has been prepared to describe a spanning tree of hits.
/// The input clusters must have the hits ordered along the branches of the
/// tree, and must start with the end point hits (front and back) being shared
/// with the neighboring clusters.  For the input clusters, all of the
/// neighbors share end point hits, and none of the neighbors share interior
/// hits.
class Cube::GrowClusters: public Cube::Algorithm {
public:
    GrowClusters();
    virtual ~GrowClusters() {}

    /// The routine that does the actual work.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// Set the maximum number of hits to consider at either end of a cluster
    /// when checking to see if it's consistent with a line.
    void SetMaxLineHits(int n) {fMaxLineHits = n;}

    /// The acceptance threshold for combining two segments.  The chi2 is for
    /// a single degree of freedom, so that can guide the approximate best
    /// value.  A value of 6.0 has a p-value of about 1.4 percent.
    void SetChi2Threshold(double c) {fChi2Threshold = c;}

private:

    /// The maximum number of hits to consider at either end of a cluster when
    /// checking to see if it's consistent with a line.
    int fMaxLineHits;

    /// Determine the maximum chi2 (1 dof) for an acceptable line fit.
    double fChi2Threshold;

    /// A minimum average charge per hit in a cluster (no ghosts).  Clusters
    /// below this threshold are not considered for merging, but are still
    /// available for later analysis.
    double fChargePerHitThreshold;

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
