#ifndef CubeHitSelection_hxx_seen
#define CubeHitSelection_hxx_seen

#include <iostream>
#include <vector>

#include <TObject.h>
#include <TNamed.h>
#include <TVector3.h>

#include "CubeHit.hxx"

namespace Cube {
    class HitSelection;
}

/// A container of HitHandle objects for the hit detector information.  This
/// is an enhanced vector that works well with ROOT.
class Cube::HitSelection : public TNamed,
                           public std::vector< Cube::Handle<Hit> > {
public:
    HitSelection(const char* name="hits",
                 const char* title="Hit Handles");
    virtual ~HitSelection();

    /// This is the usual std::vector::push_back, but enhanced to make
    /// sure that only valid hits are inserted into the HitSelection.
    virtual void push_back(const Cube::Handle<Cube::Hit>& hit);

    /// A convenience method to make sure that a hit is only added to the
    /// HitSelection once.  The AddHit method is much slower than a
    /// push_back(), so it should only be used when the hit might already be
    /// in the selection.
    virtual void AddHit(const Cube::Handle<Cube::Hit>&);

    /// A convenience method to make sure that a hit is not in a
    /// HitSelection.
    virtual void RemoveHit(const Cube::Handle<Cube::Hit>&);

    /// Print the data vector information.
    virtual void ls(Option_t* opt = "") const;

private:
    ClassDef(HitSelection,1);
};
#endif
