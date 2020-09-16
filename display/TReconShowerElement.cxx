#include "TReconShowerElement.hxx"

#include <CubeUnits.hxx>
#include <TUnitsTable.hxx>

#include <CubeHandle.hxx>
#include <CubeReconShower.hxx>
#include <CubeReconNode.hxx>
#include <CubeShowerState.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHandle.hxx>

#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoSphere.h>
#include <TGeoMatrix.h>

#include <TEveLine.h>
#include <TEveGeoShape.h>

#include <sstream>

Cube::TReconShowerElement::~TReconShowerElement() {}

Cube::TReconShowerElement::TReconShowerElement(Cube::ReconShower& shower,
                                               bool showUncertainty)
    : TEveElementList() {

    Cube::Handle<Cube::ShowerState> state = shower.GetState();
    TLorentzVector pos = state->GetPosition();
    TLorentzVector var = state->GetPositionVariance();
    TVector3 dir = state->GetDirection().Unit();
    TVector3 dvar = state->GetDirectionVariance();

    // This is used as the annotation, so it needs to be better.
    std::ostringstream title;
    title << "Shower(" << shower.GetUniqueID() << ") --"
          << " Nodes: " << shower.GetNodes().size()
          << ",  Energy Deposit: " << shower.GetEDeposit()
          << std::endl
          << "   Position:  ("
          << unit::AsString(pos.X(),std::sqrt(var.X()),"length")
          << ", "<<unit::AsString(pos.Y(),std::sqrt(var.Y()),"length")
          << ", "<<unit::AsString(pos.Z(),std::sqrt(var.Z()),"length")
          << ")";

    title << std::endl
          << "   Direction: ("
          << unit::AsString(dir.X(), dvar.X(),"direction")
          << ", " << unit::AsString(dir.Y(), dvar.Y(),"direction")
          << ", " << unit::AsString(dir.Z(), dvar.Z(),"direction")
          << ")";

    title << std::endl
          << "   Algorithm: " << shower.GetAlgorithmName()
          << " w/ goodness: " << shower.GetQuality()
          << " / " << shower.GetNDOF();

    std::cout << "shower:: " << title.str() << std::endl;

    Cube::ReconNodeContainer& nodes = shower.GetNodes();
    std::cout << "nodes::" << "Shower Nodes " << nodes.size() << std::endl;

    SetMainColor(kRed);
    std::ostringstream objName;
    objName << shower.GetName() << "(" << shower.GetUniqueID() << ")";
    SetName(objName.str().c_str());
    SetTitle(title.str().c_str());

    TEveLine* showerLine = new TEveLine(nodes.size());
    showerLine->SetName(objName.str().c_str());
    showerLine->SetTitle(title.str().c_str());
    showerLine->SetSourceObject(&shower);
    showerLine->SetLineColor(kRed);
    showerLine->SetLineStyle(1);
    showerLine->SetLineWidth(2);

    int p=0;

    // Start at the front state of the shower
    if (state) {
        TLorentzVector frontPos = state->GetPosition();
        TLorentzVector frontVar = state->GetPositionVariance();
        showerLine->SetPoint(p++, frontPos.X(), frontPos.Y(), frontPos.Z());
        std::cout << "nodes:: "
                  << "Front:"
                  << unit::AsString(frontPos.X(),
                                    std::sqrt(frontVar.X()),"length")
                  <<", "<<unit::AsString(frontPos.Y(),
                                         std::sqrt(frontVar.Y()),"length")
                  <<", "<<unit::AsString(frontPos.Z(),
                                         std::sqrt(frontVar.Z()),"length")
                  << std::endl;;
        std::cout << "nodes: "
                  << "Front Dir: "
                  << unit::AsString(state->GetDirection())
                  << std::endl;
    }

    // Draw the line down the center of the shower.
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        Cube::Handle<Cube::ShowerState> nodeState = (*n)->GetState();
        if (!nodeState) continue;
        TLorentzVector nodePos = nodeState->GetPosition();
        showerLine->SetPoint(p++, nodePos.X(), nodePos.Y(), nodePos.Z());
    }

    AddElement(showerLine);

    // Draw spheres at the nodes.
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        Cube::Handle<Cube::ShowerState> nodeState = (*n)->GetState();
        if (!nodeState) continue;
        TLorentzVector nodePos = nodeState->GetPosition();
        double nodeWidth = nodeState->GetCone();
        TEveGeoShape *nodeShape = new TEveGeoShape("showerNode");
        nodeShape->SetName(objName.str().c_str());
        nodeShape->SetTitle(title.str().c_str());
        nodeShape->SetMainColor(kRed);
        // Set the translation
        TGeoTranslation trans(nodeState->GetPosition().X(),
                              nodeState->GetPosition().Y(),
                              nodeState->GetPosition().Z());
        nodeShape->SetTransMatrix(trans);
        TGeoManager* saveGeom = gGeoManager;
        gGeoManager = nodeShape->GetGeoMangeur();
        TGeoShape* geoShape = new TGeoSphere(0.0, nodeWidth);
        nodeShape->SetShape(geoShape);
        gGeoManager = saveGeom;
        AddElement(nodeShape);
    }
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
