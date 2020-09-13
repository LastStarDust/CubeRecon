#ifndef THitChangeHandler_hxx_seen
#define THitChangeHandler_hxx_seen

#include "TVEventChangeHandler.hxx"

namespace Cube {
    class THitChangeHandler;
};

class TEveElementList;

/// Handle drawing the hits in the event.
class Cube::THitChangeHandler: public TVEventChangeHandler {
public:
    THitChangeHandler();
    ~THitChangeHandler();

    /// Draw the hits into the current scene.
    virtual void Apply();

private:

    /// The hits that are drawn.
    TEveElementList* fHitList;

};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
