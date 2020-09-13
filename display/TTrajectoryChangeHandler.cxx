#include "TTrajectoryChangeHandler.hxx"
#include "TEventDisplay.hxx"
#include "TGUIManager.hxx"

#include <TEventManager.hxx>

#include <TGeoManager.h>
#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>

#include <sstream>
#include <iostream>

Cube::TTrajectoryChangeHandler::TTrajectoryChangeHandler() {
    fTrajectoryList = new TEveElementList("g4Trajectories",
                                          "Geant4 Trajectories");
    fTrajectoryList->SetMainColor(kRed);
    fTrajectoryList->SetMainAlpha(1.0);
    gEve->AddElement(fTrajectoryList);
}

Cube::TTrajectoryChangeHandler::~TTrajectoryChangeHandler() {
}

void Cube::TTrajectoryChangeHandler::Apply() {

    fTrajectoryList->DestroyElements();

    if (!Cube::TEventDisplay::Get().GUI().GetShowTrajectoriesButton()->IsOn()) {
        return;
    }

    std::cout <<"Handle the trajectories" << std::endl;

    for (TG4TrajectoryContainer::iterator traj
             = gCubeSimEvent->Trajectories.begin();
         traj != gCubeSimEvent->Trajectories.end();
         ++traj) {
        const TG4Trajectory::TrajectoryPoints& points = traj->Points;

        std::ostringstream label;
        label << traj->GetName()
              << " (" << traj->GetInitialMomentum().E() << " MeV)";

        bool charged = true;

        TEveLine *track = new TEveLine();
        track->SetName("trajectory");
        track->SetTitle(label.str().c_str());
        if (charged) {
            track->SetLineColor(kRed);
            track->SetLineStyle(3);
        }
        else {
            track->SetLineColor(kRed+4);
            track->SetLineStyle(4);
        }

        for (std::size_t p = 0; p < points.size(); ++p) {
            track->SetPoint(p,
                            points[p].GetPosition().X(),
                            points[p].GetPosition().Y(),
                            points[p].GetPosition().Z());
        }
        fTrajectoryList->AddElement(track);
    }
}
