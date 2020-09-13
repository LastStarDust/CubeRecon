#include "CubeHits3D.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"
#include "CubeUnits.hxx"
#include "CubeInfo.hxx"
#include "CubeShareCharge.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>
#include <CubeReconCluster.hxx>

#include <algorithm>
#include <memory>
#include <set>
#include <cmath>

namespace {
    struct hitCompareHitZ {
        bool operator () (const Cube::Handle<Cube::Hit>& lhs,
                          const Cube::Handle<Cube::Hit>& rhs) {
            return lhs->GetPosition().Z() < rhs->GetPosition().Z();
        }
    };

};

Cube::Hits3D::Hits3D()
    : Cube::Algorithm("Cube::Hits3D", "Build 3D Cube Hits from Fibers") {
    fShareCharge = 1;
    fConserveChargeSum = 1;
    fLightSpeed = 200.0*unit::mm/unit::ns;
}

Cube::Hits3D::~Hits3D() { }

// Sort the fiber times, and find the average time for hits in a narrow
// window of time.  This works so that if one fiber is "out of time" for
// this hit (usually because the first photon came from a different cube),
// we still get the average for the fibers that are in coincidence.  If
// none of the fibers are in coincidence, this returns the average of all
// three fibers with a large uncertainty.
std::pair<double,double> Cube::Hits3D::HitTime(FiberTQ& fiberTQ) const {
    const double hitRes = 1.25*unit::ns/std::sqrt(12.0);

    // Find the best window for the hit time.
    const double timeWindow = 5.0*unit::ns;
    std::sort(fiberTQ.begin(),fiberTQ.end());
    FiberTQ::iterator begin = fiberTQ.begin();
    FiberTQ::iterator end = fiberTQ.end();
    if (fiberTQ.back().first-fiberTQ.front().first > timeWindow) {
        // The hit time spread is large, so find the ones in a window.
        FiberTQ::iterator begin = fiberTQ.begin();
        FiberTQ::iterator end = begin + 1;
        double bestCharge = begin->second;
        for (FiberTQ::iterator t = fiberTQ.begin(); t != fiberTQ.end(); ++t) {
            double theCharge = 0.0;
            for (FiberTQ::iterator tt = t; tt != fiberTQ.end(); ++tt) {
                if (tt->first - t->first > timeWindow) break;
                theCharge += tt->second;
                if (theCharge < bestCharge) continue;
                bestCharge = theCharge;
                begin = t;
                end = tt + 1;
            }
        }
    }

    // None of the hits are close together and time time is being take from
    // the hit with the most charge.  In this case, the uncertainty is the
    // total spread in the hit time.
    if (end - begin < 2) {
        double dt = fiberTQ.back().first - fiberTQ.front().first;
        if (dt < hitRes) dt = hitRes;
        return std::pair<double,double>(begin->first,dt);
    }

    // Now take the charge weighted average and RMS.
    double tSum = 0.0;
    double ttSum = 0.0;
    double qSum = 0.0;
    for (FiberTQ::iterator t = begin; t != end; ++t) {
        tSum += t->second*t->first;
        ttSum += t->second*t->first*t->first;
        qSum += t->second;
    }

    if (qSum > 0.0) {
        double tHit = tSum/qSum;
        ttSum = ttSum/qSum;
        double tRMS = ttSum - tHit*tHit;
        if (tRMS > 0.0) tRMS = std::sqrt(tRMS);
        else tRMS = 0.0;
        tRMS = std::max(tRMS,hitRes/std::sqrt(1.0*(end-begin)));
        return std::pair<double,double>(tHit,tRMS);
    }

    CUBE_ERROR << "Invalid hit time" << std::endl;
    throw std::runtime_error("Invalid cube hit");
    return std::pair<double,double>();
}

bool Cube::Hits3D::MakeHit(Cube::HitSelection& writableHits,
                             const Cube::Handle<Cube::Hit>& hit1,
                             const Cube::Handle<Cube::Hit>& hit2,
                             const Cube::Handle<Cube::Hit>& hit3) const {

    // Figure out the identifier for the new hit.
    int cubeId = -1;
    if (hit3) {
        cubeId = Cube::Info::Combine3DST(hit1->GetIdentifier(),
                                         hit2->GetIdentifier(),
                                         hit3->GetIdentifier());
    }
    else {
        CUBE_ERROR << "Missing a hit needed to build a 3D Hit" << std::endl;
        return false;
    }

    // Figure out the position for the new hit.
    TVector3 hitPos = hit1->GetPosition();
    TVector3 hitSize = hit1->GetSize();
    if (Cube::Info::CubeNumber(hit1->GetIdentifier()) < 0) {
        hitPos.SetX(hit2->GetPosition().X());
        hitSize.SetX(hit2->GetSize().X());
    }
    else if (Cube::Info::CubeBar(hit1->GetIdentifier()) < 0) {
        hitPos.SetY(hit2->GetPosition().Y());
        hitSize.SetY(hit2->GetSize().Y());
    }
    else if (Cube::Info::CubePlane(hit1->GetIdentifier()) < 0) {
        hitPos.SetZ(hit2->GetPosition().Z());
        hitSize.SetZ(hit2->GetSize().Z());
    }

    // Set the charge to the total charge.
    double qHit = hit1->GetCharge();
    qHit += hit2->GetCharge();
    qHit += hit3->GetCharge();

    FiberTQ fiberTQ;

    // Get the times for the hit.
    double dHit1 = (hitPos - hit1->GetPosition()).Mag();
    double dHit2 = (hitPos - hit2->GetPosition()).Mag();
    double dHit3 = (hitPos - hit3->GetPosition()).Mag();

    double tHit1 = hit1->GetTime() - dHit1/fLightSpeed;
    double tHit2 = hit2->GetTime() - dHit2/fLightSpeed;
    double tHit3 = hit3->GetTime() - dHit3/fLightSpeed;

    fiberTQ.push_back(std::make_pair(tHit1,hit1->GetCharge()));
    fiberTQ.push_back(std::make_pair(tHit2,hit2->GetCharge()));
    fiberTQ.push_back(std::make_pair(tHit3,hit3->GetCharge()));

    double dtHit = std::abs(tHit1-tHit2);
    dtHit = std::max(dtHit,std::abs(tHit1-tHit3));
    dtHit = std::max(dtHit,std::abs(tHit2-tHit3));

    // Check that the hits are "at the same time".  This is the largest
    // possible different between simultaneous cube times.  It's determined by
    // geometry.
    if (dtHit > 20*unit::ns) return false;

    // Protect against future special cases.
    if (fiberTQ.empty()) return false;

    // Create the new writable hit.
    Cube::Handle<Cube::WritableHit> hit(new Cube::WritableHit);
    hit->SetIdentifier(cubeId);
    hit->SetPosition(hitPos);
    hit->SetSize(hitSize);

    hit->SetCharge(qHit);
    hit->SetChargeUncertainty(std::sqrt(qHit));

    std::pair<double,double> timePair = HitTime(fiberTQ);
    hit->SetTime(timePair.first);
    hit->SetTimeUncertainty(timePair.second);

    // Set the position uncertainty with charge weighting so that average
    // positions are charge weighted.  This is needed so that ReconClusters
    // are correctly genererated.  The actual uncertainty for a single hit is
    // still the size of a cube.
    double dx = std::sqrt(1.0*unit::cm*unit::cm/12.0/qHit);
    TVector3 unc(dx,dx,dx);
    hit->SetUncertainty(unc);

    // Add the first hit.  There MUST be a first hit!
    hit->AddHit(hit1);
    for (int i=0; i<hit1->GetContributorCount(); ++i) {
        hit->AddContributor(hit1->GetContributor(i));
    }

    // Add the second hit.  There MUST be a second hit!
    hit->AddHit(hit2);
    for (int i=0; i<hit2->GetContributorCount(); ++i) {
        hit->AddContributor(hit2->GetContributor(i));
    }

    // Add the third hit.  There SHOULD be a third hit, but it's not required.
    hit->AddHit(hit3);
    for (int i=0; i<hit3->GetContributorCount(); ++i) {
        hit->AddContributor(hit3->GetContributor(i));
    }

    writableHits.push_back(hit);

    return true;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::Hits3D::Process(const Cube::AlgorithmResult& fibers,
                      const Cube::AlgorithmResult&,
                      const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Hits3D::Process" << std::endl;

    Cube::Handle<Cube::HitSelection> fiberHits = fibers.GetHitSelection();
    if (!fiberHits) {
        CUBE_ERROR << "No input hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    typedef std::set< Cube::Handle<Cube::Hit> >  HitSet;
    HitSet usedSet;
    HitSet unusedSet;

    Cube::HitSelection xzHits;
    Cube::HitSelection yzHits;
    std::map<int, Cube::HitSelection> xyHits;
    for (Cube::HitSelection::iterator h = fiberHits->begin();
         h != fiberHits->end(); ++h) {
        int plane = Cube::Info::IdentifierProjection((*h)->GetIdentifier());
        switch (plane) {
        case Cube::Info::kXZProj:
            xzHits.push_back(*h);
            break;
        case Cube::Info::kYZProj:
            yzHits.push_back(*h);
            break;
        case Cube::Info::kXYProj:
            xyHits[(*h)->GetIdentifier()].push_back(*h);
            break;
        default:
            CUBE_ERROR << "Invalid fiber plane" << std::endl;
        }
    }

    std::sort(xzHits.begin(), xzHits.end(), hitCompareHitZ());
    std::sort(yzHits.begin(), yzHits.end(), hitCompareHitZ());

    CUBE_LOG(0) << "XZ Hits: " << xzHits.size()
                << " YZ Hits: " << yzHits.size()
                << " XY Hits: " << xyHits.size()
                << std::endl;

    Cube::HitSelection writableHits;
    for (Cube::HitSelection::iterator xz=xzHits.begin();
         xz!=xzHits.end(); ++xz) {
        Cube::HitSelection::iterator yzBegin = yzHits.begin();
        Cube::HitSelection::iterator yzEnd = yzHits.end();

        for (Cube::HitSelection::iterator yz=yzBegin;
             yz!=yzEnd; ++yz) {
            if ((*yz)->GetPosition().Z()<(*xz)->GetPosition().Z()-5) continue;
            if ((*yz)->GetPosition().Z()>(*xz)->GetPosition().Z()+5) continue;
            int xyFiber = Cube::Info::Identifier3DST(
                Cube::Info::CubeNumber((*xz)->GetIdentifier()),
                Cube::Info::CubeBar((*yz)->GetIdentifier()),
                -1);
            std::map<int, Cube::HitSelection>::iterator xySel
                = xyHits.find(xyFiber);
            if (xySel == xyHits.end()) continue;
            Cube::HitSelection& xy = xySel->second;
            for (Cube::HitSelection::iterator h = xy.begin();
                 h != xy.end(); ++h) {
                if (MakeHit(writableHits,*xz,*yz,*h)) {
                    usedSet.insert(*xz);
                    usedSet.insert(*yz);
                    usedSet.insert(*h);
                }
            }
        }
    }
    std::cout << "Hits Generated " << writableHits.size() << std::endl;

    // Clear the selections to make sure the handles reset.
    xzHits.clear();
    yzHits.clear();
    xyHits.clear();

    // All of the 3 fiber 3D hits have been found, now check any unassociated
    // hits.  Unassociated hits are going into the unused list.
    for (Cube::HitSelection::iterator h = fiberHits->begin();
         h != fiberHits->end(); ++h) {
        if ((*h)->GetConstituentCount() > 0) {
            continue;
        }
        if (usedSet.find(*h) != usedSet.end()) continue;
        unusedSet.insert(*h);
    }

    CUBE_LOG(0) << "Total 3D Hits: " << writableHits.size() << std::endl;

    if (writableHits.size() < 1) {
        CUBE_ERROR << "No 3D hits" << std::endl;
        return Cube::Handle<Cube::AlgorithmResult>();
    }

    // Share the charge among the 3D hits so that the total charge in the
    // event is not overcounted.  This corrects for the attenuation in the
    // fiber.
    Cube::ShareCharge shareCharge;
    shareCharge.SetChargeConservation(fConserveChargeSum);
    if (fShareCharge == 1) {
        // Apply (almost) the same formalism as for the MaximumEntropy
        // version.  Share the charge by predicting the measurement in each
        // fiber based on the deposit in each cube.  The deposit in a cube is
        // constrained to be positive.  An approximation of the maximum
        // entropy prior is applied to break degeneracies.  The prior prefers
        // that cubes have the same charge (it's a very weak prior).
        shareCharge.OptimizeCubes(writableHits);
    }
    else if (fShareCharge == 2) {
        // Share the charge assuming that all of the fibers in a cube will
        // measure the same amount of deposit, constrained by the charge
        // present at the fiber being equal to the measured charge.  The
        // constraint on the charge at the fiber is applied using a Lagrange
        // multiplier.  The constraint at the cube is the likelihood that all
        // three measurements came from the same mean.
        shareCharge.ApplyConstraints(writableHits);
    }
    else if (fShareCharge == 3) {
        // Share the charge applying a Bayesian probability with a maximum
        // entropy prior.  The probability is based on predicting the
        // measurement in each fiber based on the deposit in each cube.  The
        // deposit in a cube is constrained to be positive.  The maximum
        // entropy prior is that all cubes should have the same charge (it's a
        // very weak prior).  This is impossibly slow for large events and has
        // precision problems.  The ApplyConstraints and OptimizeCubes
        // versions are better.
        shareCharge.MaximizeEntropy(writableHits);
    }

    // Copy the writable hits into the clustered hit selection;
    Cube::HitSelection clustered;
    for (Cube::HitSelection::iterator h = writableHits.begin();
         h != writableHits.end(); ++h) {
        Cube::Handle<Cube::WritableHit> hit = *h;
        Cube::Handle<Cube::Hit> newHit(new Cube::Hit(*hit));
        clustered.push_back(newHit);
    }

    // Build an object container with the hits clustered into a convenient
    // form. This is probably mostly used for display.
    Cube::Handle<Cube::ReconObjectContainer> finalObjects(
        new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    // Create a single cluster from the new 3D hits.
    Cube::Handle<Cube::ReconCluster> cluster3D
        = Cube::CreateCluster("hits3D", clustered.begin(),clustered.end());
    if (cluster3D) {
        finalObjects->push_back(cluster3D);
    }

    // Collect all of the hits into the used and unused hit selections.
    Cube::MakeUsed makeUsed(*fiberHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
