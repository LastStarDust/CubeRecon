#include "CubeEvent.hxx"
#include "CubeMakeHits3D.hxx"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVector.h>
#include <TH1F.h>

#include <iostream>
#include <sstream>
#include <exception>
#include <memory>

int main(int argc, char **argv) {
    std::cout << "CubeRecon: Hello World" << std::endl;
    int maxEntries = 1E+8; // Maximum to process.

    while (true) {
        int c = getopt(argc,argv,"n:");
        if (c<0) break;
        switch (c) {
        case 'n': {
            std::istringstream tmp(optarg);
            tmp >> maxEntries;
            break;
        }
        default: {
            std::cout << "Usage: " << std::endl;
            std::cout << "   "
                      << "-n <number>  : Process no more than"
                      << " <number> events."
                      << std::endl;
            exit(1);
        }
        }
    }

    if (argc <= optind) {
        throw std::runtime_error("Missing input file");
    }
    std::string inputName(argv[optind++]);
    std::cout << "Input Name " << inputName << std::endl;

    if (argc <= optind) {
        throw std::runtime_error("Missing output file");
    }
    std::string outputName(argv[optind++]);
    std::cout << "Output Name " << outputName << std::endl;

    std::unique_ptr<TFile> inputFile(new TFile(inputName.c_str(),"old"));
    if (!inputFile->IsOpen()) {
        throw std::runtime_error("File not open");
    }
    std::cout << "Input File " << inputFile->GetName() << std::endl;

    /// Attach to the input tree.
    TTree* inputTree = (TTree*) inputFile->Get("CubeEvents");
    if (!inputTree) {
        std::cout << "Missing the event tree" << std::endl;
        return 1;
    }
    Cube::Event *inputEvent = NULL;
    inputTree->SetBranchAddress("Event",&inputEvent);

    // Attach to the output tree.
    std::unique_ptr<TFile> outputFile(new TFile(outputName.c_str(),"recreate"));
    TTree *outputTree = new TTree("CubeEvents","Reconstructed Event");
    static Cube::Event *outputEvent = inputEvent;
    outputTree->Branch("Event",&outputEvent);

    TH1F* hitCharge = new TH1F("hitCharge", "Charge of the hits",
                               500, 0.0, 500.0);

    // Loop through the events.
    int totalEntries = inputTree->GetEntries();
    totalEntries = std::min(totalEntries,maxEntries);
    for (int entry = 0; entry < totalEntries; ++entry) {
        inputTree->GetEntry(entry);
        outputEvent = inputEvent;
        outputEvent->MakeCurrentEvent();
        CUBE_LOG(0) << "Process event " << outputEvent->GetRunId()
                    << "/" << outputEvent->GetEventId() << std::endl;

        std::unique_ptr<Cube::MakeHits3D> makeHits3D(new Cube::MakeHits3D);
        Cube::Handle<Cube::AlgorithmResult> hits3D
            = makeHits3D->Process(*outputEvent);
        outputEvent->AddAlgorithmResult(hits3D);
        outputEvent->AddHitSelection(hits3D->GetHitSelection());

        Cube::Handle<Cube::HitSelection> finalHits
            = outputEvent->GetHitSelection();
        if (finalHits) {
            TROOT::IndentLevel();
            CUBE_LOG(0) << "Final hits " << finalHits->size()
                        << std::endl;
            for (Cube::Handle<Cube::Hit>& h : (*finalHits)) {
                hitCharge->Fill(h->GetCharge());
            }
        }

        CUBE_LOG(0) << "Finished event " << outputEvent->GetRunId()
                    << "/" << outputEvent->GetEventId() << std::endl;
        outputTree->Fill();
    }

    outputFile->Write();

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
