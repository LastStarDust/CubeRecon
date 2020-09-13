#include "CubeHitSelection.hxx"

#include <TROOT.h>

#include <algorithm>

ClassImp(Cube::HitSelection);

Cube::HitSelection::HitSelection(const char* name,
                                 const char* title)
    : TNamed(name,title) { }

Cube::HitSelection::~HitSelection() { }

void Cube::HitSelection::push_back(const Cube::Handle<Cube::Hit>& hit) {
    if (!hit) {
        std::cout << "Attempting to add a NULL hit";
        throw std::runtime_error("Invalid NULL Hit");
    }
    std::vector< Cube::Handle<Cube::Hit> >::push_back(hit);
}

void Cube::HitSelection::AddHit(const Cube::Handle<Cube::Hit>& hit) {
    Cube::HitSelection::iterator location
        = std::find(begin(), end(), hit);
    if (location == end()) push_back(hit);
}

void Cube::HitSelection::RemoveHit(const Cube::Handle<Cube::Hit>& hit) {
    Cube::HitSelection::iterator location
        = std::find(begin(), end(), hit);
    if (location != end()) erase(location);
}

void Cube::HitSelection::ls(Option_t* opt) const {
    // Print the standard header
    std::string option(opt);
    TROOT::IndentLevel();
    std::cout << ClassName() << "(" << this << ")::"
              << GetName();
    if (option.find("title")!= std::string::npos) {
        std::cout << "/" << GetTitle();
    }
    if (GetUniqueID() > 0) std::cout << " uid: " << GetUniqueID();
    if (option.find("size") != std::string::npos && Class()) {
        std::cout << " size: " << Class()->Size() << " b";
    }
    std::cout << std::endl;
    // Handle the class specific stuff.
    if (option.find("dump") != std::string::npos
        || option.find("hits") != std::string::npos) {
        TROOT::IncreaseDirLevel();
        for (const_iterator v = begin();
             v != end();
             ++v) {
            v->ls(opt);
        }
        TROOT::DecreaseDirLevel();
    }
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
