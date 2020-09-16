#ifndef TFindResultsHandler_hxx_seen
#define TFindResultsHandler_hxx_seen

#include "TVEventChangeHandler.hxx"

namespace Cube {
    class TFindResultsHandler;
};

class TEveElementList;

/// Look through the AlgorithmResults saved in an event, and add them to the
/// GUI so they can be selected.
class Cube::TFindResultsHandler: public TVEventChangeHandler {
public:
    TFindResultsHandler();
    ~TFindResultsHandler();

    /// Draw fit information into the current scene.
    virtual void Apply();

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
