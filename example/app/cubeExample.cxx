#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHit.hxx>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <iostream>
#include <sstream>

void AnalyzeEvent(const Cube::Event& event) {
    Cube::Handle<Cube::HitSelection> hits
        = event.GetHitSelection();
    if (hits) {
        std::cout << "Selected Hits: " << hits->size() << std::endl;
    }

    static TH1F* histEventHits = NULL;
    if (!histEventHits) {
        std::cout << "Create the histogram" << std::endl;
        histEventHits = new TH1F("eventHits","Hits in the event",
                                 40,0.0,10000.0);
    }
    histEventHits->Fill(hits->size());

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event.GetObjectContainer();
    if (objects) {
        std::cout << "Reconstructed Objects: " << objects->size() << std::endl;
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
