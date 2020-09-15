#ifndef CubeClusterHits_hxx_seen
#define CubeClusterHits_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class ClusterHits;
};

/// Do a density clustering for the input hits based on the distance between
/// hits.  The cubes are assumed to have already been split up by type, so the
/// time difference between cubes is not used.  This expects the input hits to
/// be for cubes (i.e. 3D hits), and produces a reconstruction object
/// container of clusters.  The default clustering parameters are setup to
/// make sure that the cubes in the clusters are simply connected.
class Cube::ClusterHits: public Cube::Algorithm {
public:
    ClusterHits();
    virtual ~ClusterHits() {}

    /// The routine that does the actual work.
    Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// Set the neighborhood to check for neighbors.  This is the half size of
    /// the region to count cubes.  This is in units of "cubes", not distance,
    /// so neighbors are 1 cube apart, not 10 mm.
    void SetCubeNeighborhood(int n) {fCubeNeighborhood = n;}

    /// Set the number of cubes that need to be in the neighborhood to form a
    /// cluster.
    void SetCubeCount(int n) {fMinimumPoints = n;}

private:

    /// The "radius" of the box that defines the cube neighborhood.  This
    /// needs to be one or more.
    int fCubeNeighborhood;

    /// The number of neighbors needed to form a new neighborhood.
    int fMinimumPoints;
};
#endif
