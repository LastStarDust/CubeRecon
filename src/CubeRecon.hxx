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


    void SetOversizeCut(int i) {fOversizeCut = i;}

private:

    // The clustering implementation can be pretty slow, so protect against
    // really large number of hits [the current simple implementation is
    // O(N^2)!!!].  This keeps the code from hanging on large events.  It
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
