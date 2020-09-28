#include "TG4HitChangeHandler.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TGUIManager.hxx"

#include "CubeEvent.hxx"

#include <TEventManager.hxx>

#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>

#include <sstream>
#include <iostream>
#include <iomanip>

Cube::TG4HitChangeHandler::TG4HitChangeHandler() {
    fG4HitList = new TEveElementList("g4HitList","Geant4 Truth Hits");
    fG4HitList->SetMainColor(kCyan);
    fG4HitList->SetMainAlpha(1.0);
    gEve->AddElement(fG4HitList);
}

Cube::TG4HitChangeHandler::~TG4HitChangeHandler() {
}

void Cube::TG4HitChangeHandler::Apply() {

    fG4HitList->DestroyElements();

    if (!Cube::TEventDisplay::Get().GUI().GetShowG4HitsButton()->IsOn()) {
        std::cout << "G4 hits disabled" << std::endl;
        return;
    }

    double minEnergy = 0.10;
    double maxEnergy = 1.0;

    if (!Cube::gEvent) return;
    std::cout << "Handle the geant4 truth hits" << std::endl;

    TEveElementList* hitList
        = new TEveElementList( "hits","Energy deposit");
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator g4Hit
             = Cube::gEvent->G4Hits.begin();
         g4Hit != Cube::gEvent->G4Hits.end();
         ++g4Hit) {
        double energy = (*g4Hit)->GetEnergyDeposit();
        double length = ((*g4Hit)->GetStop()
                         - (*g4Hit)->GetStart()).Vect().Mag();
        double dEdX = energy;
        if (length>0.01) dEdX /= length;

        TEveLine* eveHit = new TEveLine(2);
        std::ostringstream title;
        title << "Hit(" << (*g4Hit)->GetSegmentId()
              << ":" << (*g4Hit)->GetPrimaryId()
              << ":"<<(*g4Hit)->GetPDG() << ")";
        title << std::fixed << std::setprecision(2)
              << " " << dEdX << " MeV/mm";
        title << " for " << length << " mm";

        eveHit->SetTitle(title.str().c_str());
        eveHit->SetLineWidth(1);
        eveHit->SetLineColor(TEventDisplay::Get().LogColor(
                                 dEdX,
                                 minEnergy,
                                 maxEnergy,
                                 3));
        eveHit->SetPoint(0,
                         (*g4Hit)->GetStart().X(),
                         (*g4Hit)->GetStart().Y(),
                         (*g4Hit)->GetStart().Z());
        eveHit->SetPoint(1,
                         (*g4Hit)->GetStop().X(),
                         (*g4Hit)->GetStop().Y(),
                         (*g4Hit)->GetStop().Z());
        hitList->AddElement(eveHit);
    }
    fG4HitList->AddElement(hitList);
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
