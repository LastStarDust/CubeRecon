#ifndef ToolTrueDirection_hxx_seen
#define ToolTrueDirection_hxx_seen

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

#include <TVector3.h>

namespace Cube {
    namespace Tool {
        TVector3
        ObjectTrueDirection(Cube::Event& event, Cube::ReconObject& object);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
