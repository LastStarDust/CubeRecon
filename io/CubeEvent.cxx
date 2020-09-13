#include "CubeEvent.hxx"
#include "CubeAlgorithmResult.hxx"
#include "CubeLog.hxx"

#include <TROOT.h>

ClassImp(Cube::Event)

namespace Cube {
    Event* gCurrentEvent = NULL;
};

Cube::Event* Cube::Event::CurrentEvent() {return Cube::gCurrentEvent;}

void Cube::Event::MakeCurrentEvent() const {
    Cube::gCurrentEvent = const_cast<Cube::Event*>(this);
}

Cube::Event::Event() {Initialize();}

Cube::Event::~Event() {}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
