#ifndef CubeVertexFit_hxx_seen
#define CubeVertexFit_hxx_seen

#include <CubeReconVertex.hxx>
#include <CubeReconTrack.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class VertexFit;
};

/// A class to fit a vertex.  The vertex is expected to have the tracks
/// already attached, and the fit will be the position that minimizes the
/// distances of closest approach.  The input vertex is expected to be
/// modified by the fitter so that the resulting handle will often be equal to
/// the input handle.  However, this is not guaranteed.  The resulting vertex
/// may be a different object than the input vertex.  If the fit fails, this
/// returns a NULL handle.
///
/// \code
/// Cube::VertexFit vertexFit;
/// Cube::Handle<ReconVertex> fittedVertex = vertexFit(inputVertex);
/// if (fittedVertex) std::cout << "fit was successful" << std::endl;
/// if (!fittedVertex) std::cout << "fit failed" << std::endl;
/// \endcode
///
/// If the fit fails then the returned handle will be empty.  The vertex fit
/// can also be applied using the Apply method which will be more convenient
/// when the fitter is referenced by a pointer.
///
/// \code
/// std::unique_ptr<Cube::VertexFit> vertexFit(new Cube::VertexFit);
/// Cube::Handle<ReconVertex> fittedVertex = vertexFit->Apply(inputVertex);
/// \endcode
///
/// \warning The input vertex is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guaranteed.  The resulting vertex may be a different object than the input
/// vertex.
class Cube::VertexFit {
public:
    explicit VertexFit();
    virtual ~VertexFit();

    /// A simple method to apply the vertex fit.  The input vertex is expected
    /// to be modified during the fit.  See Cube::VertexFit::Apply() for
    /// details.
    Cube::Handle<Cube::ReconVertex>
    operator ()(Cube::Handle<Cube::ReconVertex>& input) {return Apply(input);}

    /// Apply the fit a vertex.  The vertex is expected to have the tracks
    /// already attached, and the fit will be the position that minimizes the
    /// distances of closest approach.  The input vertex is expected to be
    /// modified by the fitter so that the resulting handle will often be
    /// equal to the input handle.  However, this is not guaranteed.  The
    /// resulting vertex may be a different object than the input vertex.  If
    /// the fit fails, this returns a NULL handle.
    virtual Cube::Handle<Cube::ReconVertex>
    Apply(Cube::Handle<Cube::ReconVertex>& input);

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
