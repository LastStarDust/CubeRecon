#ifndef CubeTimeSlice_hxx_seen
#define CubeTimeSlice_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class TimeSlice;
}

/// A algorithm that breaks the hits into groups by time.
///
/// The Process method expects a TAlgorithmResult with hits from the entire
/// event saved in the main THitSelection (e.g. input.GetHitSelection()).
/// This returns a TAlgorithmResult with a single hit selection, and one
/// reconstruction object container.  The hit selection is
///
///    - used -- All of the hits.  None are rejected.
///
/// The reconstruction object container is
///
///    - final -- This contains TReconCluster objects built from the hits
///         in used.
///
class Cube::TimeSlice: public Cube::Algorithm {
public:
    TimeSlice();
    virtual ~TimeSlice();

    Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    /// This sets the time gap between hits that is allowed.  If the gap is
    /// larger than this, then a new slice is started.
    void SetGap(double gap) {fGapCut = gap;}

private:

    double fGapCut;
};
#endif
