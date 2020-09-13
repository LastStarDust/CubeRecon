#include "CubeAlgorithm.hxx"

ClassImp(Cube::Algorithm);

Cube::Algorithm::Algorithm(const char* name, const char* title) :
    TNamed(name, title){
    SetVersion("none set");
}

Cube::Algorithm::~Algorithm(){ }

void Cube::Algorithm::Clear(Option_t*) {}

void Cube::Algorithm::SetVersion(const char* v) {
    fVersion = v;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::Algorithm::Process(const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&) {
    // The very most trivial Process method that can be written.  It does
    // nothing, and has a NULL result.
    return Cube::Handle<Cube::AlgorithmResult>();
}
