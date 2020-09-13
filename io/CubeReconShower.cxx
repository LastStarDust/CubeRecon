#include "CubeReconShower.hxx"
#include "CubeReconNode.hxx"

ClassImp(Cube::ReconShower);

Cube::ReconShower::ReconShower() {
    fState = new Cube::ShowerState;
    fNodes = new Cube::ReconNodeContainerImpl<Cube::ShowerState>;
}

Cube::ReconShower::ReconShower(const Cube::ReconShower& shower)
    : Cube::ReconObject(shower) {

    fNodes = new Cube::ReconNodeContainerImpl<Cube::ShowerState>;

    // Copy the nodes.  Create new nodes with ShowerState's
    Cube::ReconNodeContainer::const_iterator in;
    for (in=shower.GetNodes().begin(); in!=shower.GetNodes().end(); ++in) {
        Cube::Handle<Cube::ReconNode> node(new Cube::ReconNode);
        Cube::Handle<Cube::ReconObject> object = (*in)->GetObject();
        node->SetObject(object);
        Cube::Handle<Cube::ShowerState> tstate = (*in)->GetState();
        if (tstate){
            Cube::Handle<Cube::ReconState> pstate(
                new Cube::ShowerState(*tstate));
            node->SetState(pstate);
        }
        node->SetQuality((*in)->GetQuality());

        fNodes->push_back(node);
    }


    if (shower.GetState()) {
        Cube::Handle<Cube::ShowerState> state = shower.GetState();
        fState = new Cube::ShowerState(*state);
    }
    else {
        fState = new Cube::ShowerState;
    }
}

Cube::ReconShower::~ReconShower() {}

double Cube::ReconShower::GetEDeposit() const {
    Cube::Handle<Cube::ShowerState> state = GetState();
    if (!state) throw std::runtime_error("Shower missing state");
    return state->GetEDeposit();
}

TLorentzVector Cube::ReconShower::GetPosition() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ShowerState> state = GetState();
    if (!state) throw std::runtime_error("Shower missing state");
    return state->GetPosition();
}


TLorentzVector Cube::ReconShower::GetPositionVariance() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ShowerState> state = GetState();
    if (!state) throw std::runtime_error("Shower missing state");
    return state->GetPositionVariance();
}

bool Cube::ReconShower::IsXShower() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.X())) return false;
    return true;
}

bool Cube::ReconShower::IsYShower() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Y())) return false;
    return true;
}

bool Cube::ReconShower::IsZShower() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Z())) return false;
    return true;
}

int Cube::ReconShower::GetDimensions() const{
    TLorentzVector var = GetPositionVariance();
    int dim = 0;
    if (IsXShower()) ++dim;
    if (IsYShower()) ++dim;
    if (IsZShower()) ++dim;
    return dim;
}

TVector3 Cube::ReconShower::GetDirection() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ShowerState> state = GetState();
    if (!state) throw std::runtime_error("Shower missing state");
    return state->GetDirection();
}

double Cube::ReconShower::GetConeAngle() const {
    Cube::Handle<Cube::ShowerState> state = GetState();
    if (!state) throw std::runtime_error("Shower missing state");
    return state->GetCone();
}
