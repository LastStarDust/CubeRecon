#ifndef CubeMakeHits_hxx_seen
#define CubeMakeHits_hxx_seen

#include <vector>
#include <string>
#include <memory>

#include <CubeAlgorithm.hxx>
#include <CubeHitSelection.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class MakeHits3D;
};

/// This is a "pass zero" reconstruction to build the 3D hits in all of the
/// time slices.
class Cube::MakeHits3D: public Cube::Algorithm {
public:
    MakeHits3D();
    virtual ~MakeHits3D();

    /// The routine that does the actual work.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
