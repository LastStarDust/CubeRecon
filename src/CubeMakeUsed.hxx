#ifndef CubeMakeUsed_hxx_seen
#define CubeMakeUsed_hxx_seen

#include <CubeHitSelection.hxx>
#include <CubeReconObject.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeAlgorithm.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class MakeUsed;
}

/// A class to take an AlgorithmResult handle, and builds hit selections of
/// the hits that were "used" and "unused" by the result.  If a hit is used by
/// the "final" TReconObjectContainer it is considered to be "used".  The
/// input AlgorithmResult is modified.  This should be called as the last step
/// of the algorithm, and will make sure that the "used" hit selection is the
/// last final selection added to the result.  That makes it the default
/// selection for any subsequent algorithm.
class Cube::MakeUsed {
public:

    /// MakeUsed is constructed with a HitSelection of all hits that are
    /// available to algorithm.  This is usually called at the end of a
    /// Algorithm and "allHits" should be all of the input hits to the
    /// algorithm.
    explicit MakeUsed(const Cube::HitSelection& allHits);
    virtual ~MakeUsed();

    /// Take the input AlgorithmResult, and build a "used" and "unused" hit
    /// selection out of the hits in "allHits".  The new hit selections are
    /// then added to the algorithm result and returned.  The input algorithm
    /// is modified (and returned).
    Cube::Handle<Cube::AlgorithmResult>
    operator ()(Cube::Handle<Cube::AlgorithmResult> input);

private:
    /// All of the hits that might appear in the input algorithm.
    Cube::HitSelection fAllHits;
};

#endif
