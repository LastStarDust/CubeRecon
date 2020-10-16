#include "ToolCubeTruth.hxx"
#include "ToolG4Hits.hxx"

#include <CubeEvent.hxx>
#include <CubeHit.hxx>
#include <CubeHandle.hxx>
#include <CubeUnits.hxx>

#include <vector>

// Find the total energy that was deposited inside a cube.  (This should
// include crosstalk).
double Cube::Tool::CubeDeposit(Cube::Event& event,
                               Cube::Handle<Cube::Hit> hit) {
    std::vector<Cube::Handle<Cube::G4Hit>> hits
        = Cube::Tool::HitG4Hits(event,hit);

    double size = hit->GetSize().X();

    double deposit = 0.0;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator h = hits.begin();
         h != hits.end(); ++h) {
        TVector3 avg = 0.5*((*h)->GetStart().Vect() + (*h)->GetStop().Vect());
        double v = std::abs(hit->GetPosition().X() - avg.X());
        if (v > size) continue;
        v = std::abs(hit->GetPosition().Y() - avg.Y());
        if (v > size) continue;
        v = std::abs(hit->GetPosition().Z() - avg.Z());
        if (v > size) continue;
        deposit += (*h)->GetEnergyDeposit();
    }

    return deposit;
}

double Cube::Tool::CubeCrossTalk(Cube::Event& event,
                                 Cube::Handle<Cube::Hit> hit) {
    std::vector<Cube::Handle<Cube::G4Hit>> hits
        = Cube::Tool::HitG4Hits(event,hit);

    double size = hit->GetSize().X();

    double deposit = 0.0;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator h = hits.begin();
         h != hits.end(); ++h) {
        TVector3 avg = 0.5*((*h)->GetStart().Vect() + (*h)->GetStop().Vect());
        double x = std::abs(hit->GetPosition().X() - avg.X());
        double y = std::abs(hit->GetPosition().Y() - avg.Y());
        double z = std::abs(hit->GetPosition().Z() - avg.Z());
        if (x < size && y < size && z < size) continue;
        deposit += (*h)->GetEnergyDeposit();
    }

    return deposit;
}


// Find the time that the cube is "hit".  This is the time that should be
// reconstructed for the cube.  It's close to the earliest hit segment time
// (but very low energy segments are ignored).  This treats direct energy and
// cross talk the same way, so the hit time may come from the large energy
// deposition neighbor.  If the cube doesn't get over the threshold, then the
// time of the last g4 hit is used.
double Cube::Tool::CubeTime(Cube::Event& event,
                            Cube::Handle<Cube::Hit> hit,
                            double threshold) {
    std::vector<Cube::Handle<Cube::G4Hit>> hits
        = Cube::Tool::HitG4Hits(event,hit);

    std::vector<std::pair<double,double>> timeEnergy;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator h = hits.begin();
         h != hits.end(); ++h) {
        timeEnergy.push_back(std::make_pair((*h)->GetStart().T(),
                                            (*h)->GetEnergyDeposit()));
    }

    if (timeEnergy.empty()) return 1E+20;

    double enr = 0.0;
    std::sort(timeEnergy.begin(), timeEnergy.end());
    for (std::vector<std::pair<double,double>>::iterator
             tt = timeEnergy.begin();
         tt != timeEnergy.end(); ++tt) {
        enr += tt->second;
        if (enr > threshold) return tt->first;
    }

    return timeEnergy.back().first;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
