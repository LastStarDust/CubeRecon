#ifndef CubeHitUtilities_hxx_seen
#define CubeHitUtilities_hxx_seen
#include <CubeHitSelection.hxx>
#include <CubeReconObject.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    /// Collect all of the hits used by ReconObject objects in a
    /// reconstruction object container into a single hit selection.
    Cube::Handle<Cube::HitSelection>
    AllHitSelection(Cube::ReconObjectContainer& input);

    /// Copy all of the hits used by in a hitselection to a new hit selection.
    /// This is mostly for symmetry with the composite and simple hit
    /// selection functions.  This does make sure all of the hits in the
    /// selection are unique.
    Cube::Handle<Cube::HitSelection>
    AllHitSelection(Cube::HitSelection& input);

    /// Collect all of the composite hits used by ReconObject objects in a
    /// reconstruction object container into one hit selection.
    Cube::Handle<Cube::HitSelection>
    CompositeHitSelection(Cube::ReconObjectContainer& input);

    /// Collect all of the composite hits in a HitSelection into a new
    /// HitSelection.
    Cube::Handle<Cube::HitSelection>
    CompositeHitSelection(Cube::HitSelection& input);

    /// Collect all of the simple (also called fiber) hits used by
    /// ReconObject objects in a reconstruction object container into
    /// one hit selection.  This will unpack any composite hits to get
    /// the simple hits used to construct the composite hits.
    Cube::Handle<Cube::HitSelection>
    SimpleHitSelection(Cube::ReconObjectContainer& input);

    /// Collect all of the simple (also called fiber) hits in a
    /// HitSelection into a new THitSelection.  This will unpack any
    /// composite hits to get the simple hits used to construct the
    /// composite hit.
    Cube::Handle<Cube::HitSelection>
    SimpleHitSelection(Cube::HitSelection& input);
}
#endif
