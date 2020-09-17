#include "CubeMakeHits3D.hxx"
#include "CubeTimeSlice.hxx"
#include "CubeHits3D.hxx"
#include "CubeMakeUsed.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>

#include <sstream>
#include <iomanip>
#include <memory>

Cube::MakeHits3D::MakeHits3D()
    : Cube::Algorithm("MakeHits3D","Build 2D hits into 3D hits") {
}

Cube::MakeHits3D::~MakeHits3D() {
}

Cube::Handle<Cube::AlgorithmResult>
Cube::MakeHits3D::Process(const Cube::AlgorithmResult& input,
                          const Cube::AlgorithmResult&,
                          const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "MakeHits3D::Process" << std::endl;

    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        CUBE_ERROR << "No hits" << std::endl;
        return result;
    }
    if (inputHits->empty()) {
        CUBE_ERROR << "Hits empty" << std::endl;
        return result;
    }

    // Slice the event up by time.
    Cube::Handle<Cube::AlgorithmResult> timeSlice
        = Run<Cube::TimeSlice>(*inputHits);
    if (!timeSlice) return Cube::Handle<Cube::AlgorithmResult>();
    result->AddAlgorithmResult(timeSlice);

    // Create the container for the final objects.
    Cube::Handle<ReconObjectContainer> finalObjects(
        new Cube::ReconObjectContainer("final"));

    // Process each time slice as a separate cube reconstruction.
    Cube::Handle<Cube::ReconObjectContainer> slices
        = timeSlice->GetObjectContainer("final");

    int count = 0;
    for (Cube::ReconObjectContainer::iterator object = slices->begin();
         object != slices->end(); ++object) {
        Cube::Handle<Cube::HitSelection> hits = (*object)->GetHitSelection();
        if (hits && !hits->empty()) {
            // This is expecting fiber hits, so the first step is to combine
            // them into cube hits.  Notice that since this is usually
            // happening "out-of-band".  We need to make a local copy to allow
            // the THitSelection to be translated into a Cube::AlgorithmResult.
            Cube::HitSelection local("fibers");
            std::copy(hits->begin(), hits->end(), std::back_inserter(local));
            Cube::Handle<Cube::AlgorithmResult> hits3D
                = Run<Cube::Hits3D>(local);
            if (!hits3D) continue;
            std::ostringstream nm;
            nm << "Slice"
               << std::setfill('0')
               << std::setw(2)
               << count++;
            hits3D->SetName(nm.str().c_str());
            Cube::Handle<Cube::ReconObjectContainer> objects
                = hits3D->GetObjectContainer("final");
            if (objects) {
                std::copy(objects->begin(), objects->end(),
                          std::back_inserter(*finalObjects));
            }
            result->AddAlgorithmResult(hits3D);
        }
    }

    // Save the final objects last.
    result->AddObjectContainer(finalObjects);

    // Collect all of the hits into the used and unused hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    CUBE_LOG(0) << "All done" << std::endl;
    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
