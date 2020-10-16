#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeEvent.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolCubeTruth.hxx>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>

/// Filter through tracks, and assign them to trajectories.  Then check the
/// timing to see if it properly tags the track direction.
void AnalyzeEvent(Cube::Event& event) {

    static TH1F* histHitTiming = NULL;
    static TH1F* histAvgHitTiming = NULL;
    static TH1F* histMaxHitTiming = NULL;
    static TH1F* histLateHitTiming = NULL;
    if (!histHitTiming) {
        std::cout << "Create the histogram" << std::endl;
        histHitTiming = new TH1F("HitTiming","Time Resolution for Default",
                                 80,-20.0,20.0);
        histAvgHitTiming = new TH1F("AvgHitTiming","Time Resolution for Avg.",
                                    80,-20.0,20.0);
        histMaxHitTiming = new TH1F("MaxHitTiming",
                                    "Time Resolution for Last Hit",
                                    80,-20.0,20.0);
        histLateHitTiming = new TH1F("LateHitTiming",
                                     "Time Resolution for Avg. Near Last Hit",
                                     80,-20.0,20.0);
    }

    Cube::Handle<Cube::HitSelection> eventHits
        = event.GetHitSelection();
    if (!eventHits) return;

    for (Cube::HitSelection::iterator h = eventHits->begin();
         h != eventHits->end(); ++h) {
        if (!(*h)) continue;
        double energy = Cube::Tool::CubeDeposit(event,*h);
        double xtalk =  Cube::Tool::CubeCrossTalk(event,*h);
        double trueTime =  Cube::Tool::CubeTime(event,*h);
        histHitTiming->Fill((*h)->GetTime()-trueTime-100.0);
        double avgT = 0.0;
        double avgW = 0.0;
        double maxT = -1E+20;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - dd/200.0;
            avgT += tt*ch->GetCharge();
            avgW += ch->GetCharge();
            maxT = std::max(tt,maxT);
        }
        avgT = avgT/avgW;
        histMaxHitTiming->Fill(maxT-trueTime-100.0);
        histAvgHitTiming->Fill(avgT-trueTime-100.0);
        double lateT = 0.0;
        double lateW = 0.0;
        for (int i = 0; i<(*h)->GetConstituentCount();++i) {
            Cube::Handle<Cube::Hit> ch = (*h)->GetConstituent(i);
            double dd = ((*h)->GetPosition() - ch->GetPosition()).Mag();
            double tt = ch->GetTime() - dd/200.0;
            if (tt < maxT - 2.5) continue;
            lateT += tt*ch->GetCharge();
            lateW += ch->GetCharge();
        }
        lateT = lateT/lateW;
        histLateHitTiming->Fill(lateT-trueTime-100.0);
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
