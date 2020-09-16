#include "TFindResultsHandler.hxx"
#include "TEventDisplay.hxx"
#include "TEventManager.hxx"
#include "TGUIManager.hxx"

#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>

#include <TGeoManager.h>
#include <TGButton.h>

#include <TEveManager.h>
#include <TEveLine.h>

#include <TPRegexp.h>

#include <sstream>

Cube::TFindResultsHandler::TFindResultsHandler() { }

Cube::TFindResultsHandler::~TFindResultsHandler() { }

void Cube::TFindResultsHandler::Apply() {
    std::cout << "Find results" << std::endl;
    if (!Cube::gEvent) return;
    std::cout << "have event" << std::endl;

    TGTextEntry* defResult = Cube::TEventDisplay::Get().GUI().GetDefaultResult();
    TGListBox* resultsList = Cube::TEventDisplay::Get().GUI().GetResultsList();

    std::string defaultResult(defResult->GetText());
    TPRegexp regularExp(defResult->GetText());

    resultsList->RemoveAll();
    int id = 0;
    // Forage the results...
    std::vector<Cube::Handle<Cube::AlgorithmResult>> stack;
    std::vector<std::string> existingEntries;
    stack.push_back(Cube::Handle<Cube::AlgorithmResult>(Cube::gEvent,false));
    while (!stack.empty()) {
        Cube::Handle<Cube::AlgorithmResult> current = stack.back();
        std::cout << "stack entry " << stack.size() << " " << current->GetName() << std::endl;
        stack.pop_back();
        for (Cube::AlgorithmResult::ReconObjects::iterator o
                 = current->GetObjectContainers().begin();
             o !=  current->GetObjectContainers().end(); ++o) {
            std::string fullName((*o)->GetName());
            if (std::find(existingEntries.begin(),existingEntries.end(),
                          fullName) != existingEntries.end()) continue;
            existingEntries.push_back(fullName);
            resultsList->AddEntry(fullName.c_str(),++id);
            // Check to see if this result should be selected
            if (defaultResult.size() == 0) continue;
            if (!regularExp.Match(fullName.c_str())) continue;
            resultsList->Select(id);
        }
        // Add the daughter results to the stack
        for (Cube::AlgorithmResult::AlgorithmResults::const_iterator r
                 = current->GetResultsContainer().begin();
             r !=  current->GetResultsContainer().end(); ++r) {
            stack.push_back(*r);
        }
    }
    resultsList->Layout();
    resultsList->MapSubwindows();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
