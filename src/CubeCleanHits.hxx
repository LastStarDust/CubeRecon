#ifndef CubeCleanHits_hxx_seen
#define CubeCleanHits_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class CleanHits;
};

/// Do a density clustering for the input hits based on the distance between
/// hits.  The cubes are assumed to have already been split up by type, so the
/// time difference between cubes is not used.  This expects the input hits to
/// be for cubes (i.e. 3D hits), and produces a reconstruction object
/// container of clusters.  The default clustering parameters are setup to
/// make sure that the cubes in the clusters are simply connected.
class Cube::CleanHits: public Cube::Algorithm {
public:
    CleanHits();
    virtual ~CleanHits() {}

    /// The routine that does the actual work.
    Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    void SetGhostThreshold(double q) {fGhostHitThreshold = q;}

private:

    // Hits below this threshold are not considered for the reconstruction.
    double fGhostHitThreshold;
};
#endif
