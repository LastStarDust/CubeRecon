#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    static TH2F* histTrackTiming = NULL;
    static TH2F* histDirectionTiming = NULL;
    if (!histTrackTiming) {
        std::cout << "Create the histogram" << std::endl;
        histTrackTiming = new TH2F("trackTiming","Time versus length",
                                   30,0.0,3000.0,
                                   40,0.0,20.0);
        histDirectionTiming = new TH2F("directionTiming",
                                   "Time versus length rel. to direction",
                                   30,0.0,3000.0,
                                   40,-20.0,20.0);
    }

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return;

    std::cout << "Reconstructed Objects: " << objects->size() << std::endl;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
         o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) continue;
        int mainTraj = Cube::Tool::MainTrajectory(event,*track);
        if (mainTraj<0) continue;
        int primTraj = Cube::Tool::PrimaryId(event,mainTraj);
        std::vector<Cube::Handle<Cube::G4Hit>> g4Hits
            = Cube::Tool::ObjectG4Hits(event,*track);
        double trueT = 0.0;
        if (g4Hits.size()>0) {
            trueT = g4Hits.back()->GetStop().T()
                - g4Hits.front()->GetStart().T();
        }
        TVector3 trueD = Cube::Tool::ObjectTrueDirection(event,*track);
        Cube::Handle<Cube::G4Trajectory> traj = event.G4Trajectories[mainTraj];
        if (!traj) {
            std::cout << "oops" << std::endl;
        }
        TLorentzVector front = track->GetFront()->GetPosition();
        TLorentzVector back = track->GetBack()->GetPosition();
        double dT = back.T() - front.T();
        double dL = (back.Vect() - front.Vect()).Mag();
        double dCos = track->GetDirection()*trueD;
        std::cout << primTraj << " " << mainTraj << " " << trueT
                  << " " << dCos << " " << dT << " " << dL << std::endl;
        histTrackTiming->Fill(dL,dT);
        if (dCos > 0.0) histDirectionTiming->Fill(dL,dT);
        else histDirectionTiming->Fill(dL,-dT);
    }
}

int main(int argc, char** argv) {
    int maxEntries = 1E+8; // Maximum to process.
    int firstEntry = 0;

    while (true) {
        int c = getopt(argc,argv,"n:s:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
            break;
        }
        case 's': {
            std::istringstream tmp(optarg);
            tmp >> firstEntry;
            break;
        }
        default: {
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
                      << "-s <number>  : Skip <number> entries"
                      << std::endl
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl;
            exit(1);
        }
        }
    }

    if (argc <= optind) throw std::runtime_error("Missing input file");
    std::string inputName(argv[optind++]);
    std::cout << "Input Name " << inputName << std::endl;

    std::string outputName;
    if (argc > optind) {
        outputName = argv[optind++];
    }
    else {
        std::cout << "NO OUTPUT FILE!!!!" << std::endl;
    }

    // Attach to the input tree.
    std::unique_ptr<TFile> inputFile(new TFile(inputName.c_str(),"old"));
    if (!inputFile->IsOpen()) throw std::runtime_error("Input file not open");

    /// Attach to the input tree.
    TTree* inputTree = (TTree*) inputFile->Get("CubeEvents");
    if (!inputTree) throw std::runtime_error("Missing the event tree");
    Cube::Event *inputEvent = NULL;
    inputTree->SetBranchAddress("Event",&inputEvent);

    // Open the output file
    std::unique_ptr<TFile> outputFile;
    if (!outputName.empty()) {
        std::cout << "Open Output File: " << outputName << std::endl;
        outputFile.reset(new TFile(outputName.c_str(),"recreate"));
    }

    // Loop through the events.
    int totalEntries = inputTree->GetEntries();
    totalEntries = std::min(totalEntries,firstEntry+maxEntries);
    for (int entry = firstEntry; entry < totalEntries; ++entry) {
        inputTree->GetEntry(entry);
        std::cout << "Process event " << inputEvent->GetRunId()
                    << "/" << inputEvent->GetEventId() << std::endl;
        AnalyzeEvent(*inputEvent);
    }

    if (outputFile) {
        outputFile->Write();
        outputFile->Close();
    }

    return 0;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
