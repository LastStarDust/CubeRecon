#ifndef TFitChangeHandler_hxx_seen
#define TFitChangeHandler_hxx_seen

#include "TVEventChangeHandler.hxx"

#include <CubeReconObject.hxx>
#include <CubeReconCluster.hxx>
#include <CubeReconShower.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconVertex.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    class TFitChangeHandler;
};

class TEveElementList;

/// Handle drawing the TAlgorithmResults saved in the event.
class Cube::TFitChangeHandler: public TVEventChangeHandler {
public:

    TFitChangeHandler();
    ~TFitChangeHandler();

    /// Draw fit information into the current scene.
    virtual void Apply();

private:

    /// A method to draw a TReconCluster.
    int ShowReconCluster(TEveElementList* list,
                         const Cube::Handle<Cube::ReconCluster> obj,
                         int index,
                         bool forceUncertainty);

    /// A method to draw a TReconShower.
    int ShowReconShower(TEveElementList* list,
                        const Cube::Handle<Cube::ReconShower> obj,
                        int index);

    /// A method to draw a TReconTrack.
    int ShowReconTrack(TEveElementList* list,
                       const Cube::Handle<Cube::ReconTrack> obj,
                       int index);

    /// A method to draw a TReconVertex
    int ShowReconVertex(TEveElementList* list,
                        const Cube::Handle<Cube::ReconVertex> obj,
                        int index);

    /// A method to draw a generic TReconObject
    int ShowReconObject(TEveElementList* list,
                        const Cube::Handle<Cube::ReconObject> obj,
                        int index,
                        bool forceUncertainty);

    /// A method to draw a TReconObjectContainer
    int ShowReconObjects(TEveElementList* list,
                         const Cube::Handle<Cube::ReconObjectContainer> obj,
                         int index = 0);

    /// The reconstruction objects to draw in the event.
    TEveElementList* fFitList;

    /// The hits to draw in the event.
    TEveElementList* fHitList;

    /// A boolean to flag if hits associated with the fit should be drawn.
    bool fShowFitsHits;

    /// A boolean to flag if the objects associated with the fit should be
    /// shown.
    bool fShowFitsObjects;

    /// The new camera center.
    TVector3 fCameraCenter;

    /// A weight used to calculate the camera center.
    double fCameraWeight;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
