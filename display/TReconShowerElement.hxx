#ifndef TReconShowerElement_hxx_seen
#define TReconShowerElement_hxx_seen

#include <CubeReconShower.hxx>

#include <TEveElement.h>

namespace Cube {
    class TReconShowerElement;
};


/// A Eve Element object to represent a Cube::ReconShower.
class Cube::TReconShowerElement: public TEveElementList {
public:
    TReconShowerElement(Cube::ReconShower& shower, bool showUncertainty);
    virtual ~TReconShowerElement();
};
#endif
