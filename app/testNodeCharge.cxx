#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeTrackState.hxx>
#include <CubeReconCluster.hxx>
#include <CubeReconNode.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

static TH2F* histClusterCharge = NULL;
static TH2F* histStateDeposit = NULL;
static TH2F* histProtonContainedCharge = NULL;
static TH2F* histProtonContainedDeposit = NULL;
static TH2F* histProtonContainedDiff = NULL;
static TH1F* histProtonContainedProj = NULL;
static TH2F* histProtonExitingDeposit = NULL;
static TH2F* histMuonContainedDeposit = NULL;
static TH2F* histMuonExitingDeposit = NULL;
static TH1F* histMuonExitingProj = NULL;
static TH2F* histOtherContainedDeposit = NULL;
static TH2F* histOtherExitingDeposit = NULL;
static TH2F* histTrackDiff = NULL;

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {
    int backNodes = 20;
    if (!histClusterCharge) {
        std::cout << "Create the histogram" << std::endl;
        histClusterCharge = new TH2F("clusterCharge","Charge versus node",
                                   backNodes,0.0,1.0*backNodes,
                                   100,0.0,2500.0);
        histStateDeposit = new TH2F("stateDeposit","Deposit versus node",
                                    backNodes,0.0,1.0*backNodes,
                                    100,0.0,2500.0);
        histProtonContainedCharge
            = new TH2F("protonCharge",
                       "Charge versus node for contained protons",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histProtonContainedDeposit
            = new TH2F("protonContained",
                       "Deposit versus node for contained protons",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histProtonContainedDiff
            = new TH2F("protonContainedDiff",
                       "Fraction difference versus node for contained protons",
                       backNodes,0.0,1.0*backNodes,
                       40,-1.0,3.0);
        histProtonContainedProj
            = new TH1F("protonContainedDiffProj",
                       "Fraction difference versus node for contained protons",
                       40,-1.0,3.0);
        histProtonExitingDeposit
            = new TH2F("protonExiting",
                       "Deposit versus node for exiting protons",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histMuonContainedDeposit
            = new TH2F("muonContained",
                       "Deposit versus node for contained muons",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histMuonExitingDeposit
            = new TH2F("muonExiting",
                       "Deposit versus node for exiting muons",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histMuonExitingProj
            = new TH1F("muonExitingProj",
                       "Fraction difference versus node for exiting muons",
                       40,-1.0,3.0);
        histOtherContainedDeposit
            = new TH2F("otherContained",
                       "Deposit versus node for contained others",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histOtherExitingDeposit
            = new TH2F("otherExiting",
                       "Deposit versus node for exiting others",
                       backNodes,0.0,1.0*backNodes,
                       100,0.0,2500.0);
        histTrackDiff = new TH2F("trackDiff","Charge difference versus node",
                                 backNodes,0.0,1.0*backNodes,
                                 40,-1.0,3.0);
    }

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (!objects) return;

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
        double dCos = track->GetDirection()*trueD;
        int containment = Cube::Tool::ContainedObject(*track);
        if (dCos < 0.0) continue;
        Cube::Handle<Cube::G4Trajectory> traj = event.G4Trajectories[mainTraj];
        if (!traj) continue;
        double mainPurity = Cube::Tool::MainPurity(event,*track);
        if (mainPurity < 0.9) continue;
        Cube::ReconNodeContainer& nodes = track->GetNodes();
        int mainPDG = traj->GetPDGCode();
        if (nodes.size() < backNodes) continue;
        Cube::Handle<Cube::TrackState> last = nodes.back()->GetState();
        double stateDeposit = 0.0;
        double clusterDeposit = 0.0;
        for (int i=0; i<nodes.size(); ++i) {
            Cube::Handle<Cube::TrackState> st = nodes[i]->GetState();
            Cube::Handle<Cube::ReconCluster> cl = nodes[i]->GetObject();
            stateDeposit += st->GetEDeposit();
            clusterDeposit += cl->GetEDeposit();
        }
        for (int i = 1; i < nodes.size(); ++i) {
            int j = nodes.size() - i;
            Cube::Handle<Cube::TrackState> st = nodes[j]->GetState();
            Cube::Handle<Cube::ReconCluster> cl = nodes[j]->GetObject();
            if (!st) {
                std::cout << "Not a track" << std::endl;
                continue;
            }
            if (!cl) {
                std::cout << "Not a cluster" << std::endl;
                continue;
            }
            if (i > backNodes) continue;
            histStateDeposit->Fill(i-0.5, st->GetEDeposit());
            histClusterCharge->Fill(i-0.5, cl->GetEDeposit());
            double delta = cl->GetEDeposit()-st->GetEDeposit();
            delta /= st->GetEDeposit();
            histTrackDiff->Fill(i-0.5, delta);
            if (std::abs(mainPDG) == 2212) {
                if (containment > 5) {
                    histProtonContainedCharge->Fill(i-0.5, cl->GetEDeposit());
                    histProtonContainedDeposit->Fill(i-0.5, st->GetEDeposit());
                    histProtonContainedDiff->Fill(i-0.5, delta);
                    histProtonContainedProj->Fill(delta);
                }
                else if (containment < 1) {
                    histProtonExitingDeposit->Fill(i-0.5, st->GetEDeposit());
                }
            }
            else if (std::abs(mainPDG) == 13) {
                if (containment > 5) {
                    histMuonContainedDeposit->Fill(i-0.5, st->GetEDeposit());
                }
                else if (containment < 1) {
                    histMuonExitingDeposit->Fill(i-0.5, st->GetEDeposit());
                    histMuonExitingProj->Fill(delta);
                }
            }
            else {
                if (containment > 5) {
                    histOtherContainedDeposit->Fill(i-0.5, st->GetEDeposit());
                }
                else if (containment < 1) {
                    histOtherExitingDeposit->Fill(i-0.5, st->GetEDeposit());

                }
            }
        }
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
