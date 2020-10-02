#ifndef ToolMainTrajectory_hxx_seen
#define ToolMainTrajectory_hxx_seen

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

namespace Cube {
    namespace Tool {
        /// Find the trajectory that contributed most to the track.  The
        /// longest trajectory wins.
        int MainTrajectory(Cube::Event& event, Cube::ReconObject& object);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
