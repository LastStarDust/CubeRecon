#include "CubeMakeUsed.hxx"
#include "CubeHitUtilities.hxx"

#include <CubeAlgorithmResult.hxx>
#include <CubeHitSelection.hxx>

#include <algorithm>
#include <vector>
#include <set>

Cube::MakeUsed::MakeUsed(const Cube::HitSelection& allHits)
    : fAllHits(allHits) { }

Cube::MakeUsed::~MakeUsed() {}

Cube::Handle<Cube::AlgorithmResult>
Cube::MakeUsed::operator () (Cube::Handle<Cube::AlgorithmResult> input) {
    // Check for the existing unused hits.  Create the container if they don't
    // exist.
    Cube::Handle<Cube::HitSelection> unusedHits
        = input->GetHitSelection("unused");
    if (!unusedHits) {
        unusedHits = Cube::Handle<Cube::HitSelection>(
            new Cube::HitSelection("unused"));
        input->AddHitSelection(unusedHits);
        unusedHits = input->GetHitSelection("unused");
    }

    // Check for the existing used hits.  Create the container if they don't
    // exist.
    Cube::Handle<Cube::HitSelection> usedHits = input->GetHitSelection("used");
    if (!usedHits) {
        usedHits = Cube::Handle<Cube::HitSelection>(
            new Cube::HitSelection("used"));
        input->AddHitSelection(usedHits);
        usedHits = input->GetHitSelection("used");
    }

    // Get final objects from the algorithm result and add the hits to the
    // used selection.
    Cube::Handle<Cube::ReconObjectContainer> finalObjects
        = input->GetObjectContainer("final");
    if (finalObjects) {
        Cube::Handle<Cube::HitSelection> hits
            = Cube::AllHitSelection(*finalObjects);
        // The hits exist, so copy them into the usedHits.
        if (hits) {
            std::copy(hits->begin(), hits->end(),
                      std::back_inserter(*usedHits));
        }
    }

    // Make sure the used hits are unique.
    Cube::HitSelection::iterator end
        = std::unique(usedHits->begin(), usedHits->end());
    usedHits->erase(end, usedHits->end());

    // Build a selection that has all of the composite hits, as well as any
    // simple hits that were used to build the composites, plus an simple hits
    // that are directly in the used hits.
    Cube::HitSelection everyLastHit;
    std::copy(usedHits->begin(), usedHits->end(),
              std::back_inserter(everyLastHit));

    Cube::Handle<Cube::HitSelection> allTheSimpleHits
        = Cube::SimpleHitSelection(everyLastHit);
    std::copy(allTheSimpleHits->begin(), allTheSimpleHits->end(),
              std::back_inserter(everyLastHit));

    // Make sure the used hits are unique and sorted.
    end = std::unique(everyLastHit.begin(), everyLastHit.end());
    everyLastHit.erase(end, everyLastHit.end());
    std::sort(everyLastHit.begin(),everyLastHit.end());

    // Figure out which inputs hits are not used.
    if (!fAllHits.empty()) {
        // The input hits possibly a include mix of simple and composite hits.
        for (Cube::HitSelection::iterator hit = fAllHits.begin();
             hit != fAllHits.end(); ++hit) {
            if (std::binary_search(
                    everyLastHit.begin(),everyLastHit.end(),*hit)) {
                continue;
            }
            unusedHits->push_back(*hit);
        }
    }

    // Make sure the unused hits are unique.  They should be, but check
    // anyway.
    end = std::unique(unusedHits->begin(), unusedHits->end());
    unusedHits->erase(end, unusedHits->end());

    CUBE_LOG(1) << input->GetName()
                << " -- Used Hits: " << usedHits->size()
                << "   Unused Hits: " << unusedHits->size()
                << std::endl;;

    return input;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
