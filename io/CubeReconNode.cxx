
#include <TROOT.h>
#include "CubeReconNode.hxx"

ClassImp(Cube::ReconNode);

Cube::ReconNode::ReconNode()
    : fQuality(0.) { }

Cube::ReconNode::~ReconNode() {}

void Cube::ReconNode::ls(Option_t *opt) const {
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
    // Do the class specific stuff.
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << "Quality: " << GetQuality() << std::endl;
    if (fState) {
        fState->ls(opt);
    }
    else {
        TROOT::IndentLevel();
        std::cout << "State Information is missing" << std::endl;
    }
    if (fObject) {
        if (option.find("dump") != std::string::npos) {
            fObject->ls(opt);
        }
        else if (option.find("recon") != std::string::npos) {
            // Dump the top level object, but not any children.
            fObject->ls("");
        }
        else {
            TROOT::IndentLevel();
            std::cout << "Object Information not shown" << std::endl;
        }
    }
    else {
        TROOT::IndentLevel();
        std::cout << "Object Information is missing" << std::endl;
    }
    TROOT::DecreaseDirLevel();
}

ClassImp(Cube::ReconNodeContainer);

Cube::ReconNodeContainer::ReconNodeContainer() {}

Cube::ReconNodeContainer::~ReconNodeContainer() {}

void Cube::ReconNodeContainer::ls(Option_t *opt) const {
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
    // Do the class specific stuff.
    TROOT::IncreaseDirLevel();
    for (Cube::ReconNodeContainer::const_iterator n = begin();
         n != end();
         ++n) {
        (*n)->ls(opt);
    }
    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
