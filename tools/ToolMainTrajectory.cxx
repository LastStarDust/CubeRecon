#include <ToolMainTrajectory.hxx>
#include <ToolG4Hits.hxx>

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

#include <map>
#include <vector>

// Find the trajectory that contributed most to the track.  The longest
// trajectory wins.
int Cube::Tool::MainTrajectory(Cube::Event& event,
                   Cube::ReconObject& object) {
    std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
        = Cube::Tool::ObjectG4Hits(event,object);
    std::map<int,double> trajMap;
    for(std::vector<Cube::Handle<Cube::G4Hit>>::iterator
            g = g4Hits.begin(); g != g4Hits.end(); ++g) {
        trajMap[(*g)->GetPrimaryId()]
            += ((*g)->GetStart().Vect() - (*g)->GetStop().Vect()).Mag();
    }
    int maxTraj = -1;
    double maxLen = 0.0;
    for(std::map<int,double>::iterator t = trajMap.begin();
        t != trajMap.end(); ++t) {
        if (maxLen > t->second) continue;
        maxTraj = t->first;
        maxLen = t->second;
    }
    return maxTraj;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
