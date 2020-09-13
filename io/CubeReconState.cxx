#include <iostream>
#include <iomanip>
#include <cmath>

#include <TROOT.h>

#include "CubeLog.hxx"
#include "CubeReconState.hxx"
#include "CubeUnits.hxx"
#include "TUnitsTable.hxx"

ClassImp(Cube::ReconState);

Cube::ReconState::ReconState() { }

Cube::ReconState::ReconState(const ReconState& state)
    : TObject(state), fValues(state.fValues),
      fFieldNames(state.fFieldNames) { }

Cube::ReconState::~ReconState() { }

std::string Cube::ReconState::GetStateFields(void) const {
    // Construct a type name out of the field names.  This in turn is used by
    // the CorrValues class to construct a type hash which is used to make
    // sure that operations are done on compatible CorrValues objects.
    std::string typeName;
    for (std::vector<std::string>::const_iterator n = fFieldNames.begin();
         n != fFieldNames.end();
         ++n) {
        typeName += *n;
        typeName += " ";
    };
    return typeName;
}

// Build the internal state vector.
void Cube::ReconState::Init() {
    fValues.ResizeTo(fFieldNames.size());
    fValues.SetType(GetStateFields().c_str());
}

int Cube::ReconState::GetDimensions() const {
    return fValues.GetDimensions();
}

double Cube::ReconState::GetValue(int i) const {
    return fValues.GetValue(i);
}

void Cube::ReconState::SetValue(int i, double val) {
    return fValues.SetValue(i, val);
}

double Cube::ReconState::GetCovarianceValue(int i, int j) const {
    return fValues.GetCovarianceValue(i,j);
}

void Cube::ReconState::SetCovarianceValue(int i, int j, double val) {
    fValues.SetCovarianceValue(i,j,val);
}

void Cube::ReconState::SetFree(int i) {
    fValues.SetFree(i);
}

bool Cube::ReconState::IsFree(int i) const {
    return fValues.IsFree(i);
}

bool Cube::ReconState::IsFree(double v) const {
    return fValues.IsFree(v);
}

void Cube::ReconState::SetFixed(int i) {
    fValues.SetFixed(i);
}

bool Cube::ReconState::IsFixed(int i) const {
    return fValues.IsFixed(i);
}

bool Cube::ReconState::IsFixed(double v) const {
    return fValues.IsFixed(v);
}

void Cube::ReconState::Validate() {
    fValues.Validate(true);
}

Cube::CorrValues Cube::ReconState::ProjectState(
    const Cube::Handle<Cube::ReconState>& state) {
    return state->fValues;
}

/// Print the object information.
void Cube::ReconState::ls(Option_t* opt) const {
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
    // Print the class specific stuff.
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();
    std::cout << GetStateFields() << std::endl;
    TROOT::IncreaseDirLevel();
    std::ios::fmtflags save = std::cout.flags();
    for (int i = 0; i<GetDimensions(); ++i) {
        if (IsFree(i)) continue;
        TROOT::IndentLevel();
        std::cout << "  " << std::setw(6) << fFieldNames[i];
        std::cout << ":: "
                  << unit::AsString(GetValue(i), GetCovarianceValue(i,i));
        for (int j=0 ; j<i; ++j) {
            if (IsFree(j)) continue;
            if (IsFixed(j)) {
                std::cout << "  fixed";
                continue;
            }
            double c = GetCovarianceValue(i,j);
            c /= std::sqrt(GetCovarianceValue(i,i));
            c /= std::sqrt(GetCovarianceValue(j,j));
            if (std::abs(c) < 0.01) {
                std::cout << "   negl";
                continue;
            }
            std::cout << " "
                      << std::setw(6) << std::setprecision(2)
                      << c;
        }
        std::cout << std::endl;
    }
    std::cout.flags(save);
    TROOT::DecreaseDirLevel();
    TROOT::DecreaseDirLevel();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
