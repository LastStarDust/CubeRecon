#ifndef TEventManager_hxx_seen
#include <CubeEvent.hxx>

#include <TTree.h>

namespace Cube {
    /// The tree containing the event.
    extern TTree* gCubeReconTree;

    /// The event being displayed.
    extern Cube::Event *gEvent;

    /// The top algorithm in the event to display.
    extern Cube::AlgorithmResult *gAlgorithm;


    extern int gCubeReconEntryNumber;
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
