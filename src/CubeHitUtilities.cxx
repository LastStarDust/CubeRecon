#include "CubeHitUtilities.hxx"

#include <CubeReconObject.hxx>
#include <CubeHitSelection.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>

#include <set>

/// WARNING: There is a lot of duplicated code here, and it could probably be
/// implemented with some clever templating.  BUT, it's not worth changing
/// since it works.  HOWEVER, if a bug is found, you need to apply the fix to
/// each of the specialized routines.

Cube::Handle<Cube::HitSelection>
Cube::AllHitSelection(Cube::ReconObjectContainer& input) {

    // Collect every hit found into a set.
    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::ReconObjectContainer::iterator o = input.begin();
         o != input.end();
         ++o) {
        // Add the hits for the object.
        Cube::Handle<Cube::HitSelection> objHits = (*o)->GetHitSelection();
        if (objHits) {
            for (Cube::HitSelection::iterator hit = objHits->begin();
                 hit != objHits->end();
                 ++hit) {
                hitSet.insert(*hit);
            }
        }

        // Add hits for the constituents.
        Cube::Handle<Cube::ReconObjectContainer> parts
            = (*o)->GetConstituents();
        if (parts) {
            objHits = Cube::AllHitSelection(*parts);
            if (objHits) {
                for (Cube::HitSelection::iterator hit = objHits->begin();
                     hit != objHits->end();
                     ++hit) {
                    hitSet.insert(*hit);
                }
            }
        }
    }

    // Copy the set of hits into a hit selection
    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("all"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

Cube::Handle<Cube::HitSelection>
Cube::AllHitSelection(Cube::HitSelection& input) {
    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::HitSelection::iterator hit = input.begin();
         hit != input.end();
         ++hit) {
        hitSet.insert(*hit);
    }

    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("all"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

Cube::Handle<Cube::HitSelection>
Cube::CompositeHitSelection(Cube::ReconObjectContainer& input) {

    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::ReconObjectContainer::iterator o = input.begin();
         o != input.end();
         ++o) {
        // Add the hits for the object.
        Cube::Handle<Cube::HitSelection> objHits = (*o)->GetHitSelection();
        if (objHits) {
            for (Cube::HitSelection::iterator hit = objHits->begin();
                 hit != objHits->end();
                 ++hit) {
                if ((*hit)->GetConstituentCount() < 1) continue;
                hitSet.insert(*hit);
            }
        }

        // Add hits for the constituents.
        Cube::Handle<Cube::ReconObjectContainer> parts
            = (*o)->GetConstituents();
        if (parts) {
            objHits = Cube::AllHitSelection(*parts);
            if (objHits) {
                for (Cube::HitSelection::iterator hit = objHits->begin();
                     hit != objHits->end();
                     ++hit) {
                    if ((*hit)->GetConstituentCount() < 1) continue;
                    hitSet.insert(*hit);
                }
            }
        }
    }

    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("composite"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

Cube::Handle<Cube::HitSelection>
Cube::CompositeHitSelection(Cube::HitSelection& input) {
    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::HitSelection::iterator hit = input.begin();
         hit != input.end();
         ++hit) {
        if ((*hit)->GetConstituentCount() < 1) continue;
        hitSet.insert(*hit);
    }

    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("composite"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

Cube::Handle<Cube::HitSelection>
Cube::SimpleHitSelection(Cube::ReconObjectContainer& input) {
    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::ReconObjectContainer::iterator o = input.begin();
         o != input.end();
         ++o) {
        // Add the hits for the object.
        Cube::Handle<Cube::HitSelection> objHits = (*o)->GetHitSelection();
        if (objHits) {
            Cube::Handle<Cube::HitSelection> simpleHits =
                Cube::SimpleHitSelection(*(*o)->GetHitSelection());
            for (Cube::HitSelection::iterator hit = simpleHits->begin();
                 hit != simpleHits->end();
                 ++hit) {
                hitSet.insert(*hit);
            }
        }

        // Add hits for the constituents.
        Cube::Handle<Cube::ReconObjectContainer> parts
            = (*o)->GetConstituents();
        if (parts) {
            objHits = Cube::SimpleHitSelection(*parts);
            if (objHits) {
                for (Cube::HitSelection::iterator hit = objHits->begin();
                     hit != objHits->end();
                     ++hit) {
                    hitSet.insert(*hit);
                }
            }
        }
    }

    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("simple"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

Cube::Handle<Cube::HitSelection>
Cube::SimpleHitSelection(Cube::HitSelection& input) {
    std::set<Cube::Handle<Cube::Hit>> hitSet;
    for (Cube::HitSelection::iterator hit = input.begin();
         hit != input.end();
         ++hit) {
        if ((*hit)->GetConstituentCount() < 1) {
            hitSet.insert(*hit);
            continue;
        }
        // The hit is composite, so split it up.
        Cube::HitSelection compositeHits;
        for (int i=0; i < (*hit)->GetConstituentCount(); ++i) {
            compositeHits.push_back((*hit)->GetConstituent(i));
        }
        Cube::Handle<Cube::HitSelection> newHits =
            SimpleHitSelection(compositeHits);
        for (Cube::HitSelection::iterator simple = newHits->begin();
             simple != newHits->end();
             ++simple) {
            hitSet.insert(*simple);
        }
    }

    Cube::Handle<Cube::HitSelection> hits(new Cube::HitSelection("simple"));
    for (std::set<Cube::Handle<Cube::Hit>>::iterator h = hitSet.begin();
         h != hitSet.end(); ++h) {
        hits->push_back(*h);
    }

    return hits;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
