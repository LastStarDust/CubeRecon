#ifndef CubeRecon_hxx_seen
#define CubeRecon_hxx_seen

#include <CubeAlgorithm.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class Recon;
};

/// The main Reconstruction Class.  This is the top level algorithm that runs
/// the rest of the reconstruction.  It expects the 3D hits have already been
/// constructed and are provided as input.
class Cube::Recon: public Cube::Algorithm {
public:
    Recon();
    virtual ~Recon() {}

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
