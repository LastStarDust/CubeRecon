#ifndef ToolPrimaryId_hxx_seen
#define ToolPrimaryId_hxx_seen

#include <CubeEvent.hxx>
namespace Cube {
    namespace Tool {
        /// Find the primary trajectory id (this is the trajectory that was
        /// started by a GEANT4 primary particle).  The trajectories are a
        /// tree of parents with their children.  This looks at the ancestors
        /// to find the originating parent trajectory.
        int PrimaryId(Cube::Event& event, int trajId);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
