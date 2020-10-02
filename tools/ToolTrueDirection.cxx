#include "ToolTrueDirection.hxx"
#include "ToolMainTrajectory.hxx"
#include "ToolG4Hits.hxx"

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>
#include <CubeHandle.hxx>
#include <CubeG4Hit.hxx>

#include <TVector3.h>

#include <vector>

TVector3
Cube::Tool::ObjectTrueDirection(Cube::Event& event, Cube::ReconObject& object) {
    int mainTraj = Cube::Tool::MainTrajectory(event,object);
    std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
        = Cube::Tool::ObjectG4Hits(event,object);
    if (g4Hits.empty()) return TVector3();
    double maxLength = 0;
    TLorentzVector start;
    TLorentzVector stop;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator g1
             = g4Hits.begin(); g1 != g4Hits.end(); ++g1) {
        if ((*g1)->GetPrimaryId() != mainTraj) continue;
        for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator g2
                 = g4Hits.begin(); g2 != g4Hits.end(); ++g2) {
            if ((*g2)->GetPrimaryId() != mainTraj) continue;
            double length
                = ((*g2)->GetStop().Vect()-(*g1)->GetStart().Vect()).Mag();
            if (length < maxLength) continue;
            maxLength = length;
            start = (*g1)->GetStart();
            stop = (*g2)->GetStop();
        }
    }
    double dT = stop.T() - start.T();
    if (dT > 0) {
        return (stop.Vect()-start.Vect()).Unit();
    }
    return (start.Vect()-stop.Vect()).Unit();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
