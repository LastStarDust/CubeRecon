#include "ToolPrimaryId.hxx"

#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>

// Find the primary trajectory id (this is the trajectory that was started
// by a GEANT4 primary particle).
int Cube::Tool::PrimaryId(Cube::Event& event, int trajId) {
    Cube::Handle<Cube::G4Trajectory> traj = event.G4Trajectories[trajId];
    if (!traj) return -1;
    if (traj->GetParentId() < 0) return traj->GetTrackId();
    return PrimaryId(event,traj->GetParentId());
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
