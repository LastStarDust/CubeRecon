#ifndef TReconClusterElement_hxx_seen
#define TReconClusterElement_hxx_seen

#include <CubeReconCluster.hxx>

#include <TEveElement.h>

namespace Cube {
    class TReconClusterElement;
};


/// A Eve Element object to represent a Cube::ReconCluster.
class Cube::TReconClusterElement: public TEveElementList {
public:
    TReconClusterElement(Cube::ReconCluster& cluster, bool showUncertainty);
    virtual ~TReconClusterElement();

    bool Valid() {return fValid;}

private:
    bool fValid;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
