#ifndef TShowHits_hxx_seen
#define TShowHits_hxx_seen

#include <CubeHitSelection.hxx>
#include <CubeUnits.hxx>

namespace Cube {
    class TShowHits;
};

class TEveElementList;

/// Add a set of hits from a hit selection to the provided elements list.  If
/// there is a problem, this will return false.
class Cube::TShowHits {
public:

    /// Construct an object to show the Fiber hits.
    explicit TShowHits();

    /// Show the hits in the selection.  The element list is mutated by adding
    /// elements that will actually show the hit (nominally, this adds a box
    /// set).
    bool operator () (TEveElementList* elements,
                      const Cube::HitSelection& hits);

    /// Only hits with more than this charge are shown.
    void SetThreshold(double threshold) {fThreshold = threshold;}

private:

    double fThreshold;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
