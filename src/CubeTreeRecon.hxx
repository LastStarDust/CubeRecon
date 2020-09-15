#ifndef TSFGTreeRecon_hxx_seen
#define TSFGTreeRecon_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class TreeRecon;
};

/// This takes a result with a collection of 3D hits, and then breaks them
/// into clusters that are split based on the branches of a minimal spanning
/// tree (using the distance between hits).  The input hits are assumed to all
/// be in the same time slice.
class Cube::TreeRecon: public Cube::Algorithm {
public:
    TreeRecon();

    /// The routine that does the actual work.
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty);

    virtual ~TreeRecon() {}

private:
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
