#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <TROOT.h>

#include "CubeReconCluster.hxx"
#include "CubeReconNode.hxx"

ClassImp(Cube::ReconCluster);

Cube::ReconCluster::ReconCluster()
    : fMoments(3), fTemporariesInitialized(false) {
    fState = new Cube::ClusterState;
    fNodes = new Cube::ReconNodeContainerImpl<Cube::ClusterState>;
}

Cube::ReconCluster::ReconCluster(const Cube::ReconCluster& cluster)
    : Cube::ReconObject(cluster), fMoments(3), fTemporariesInitialized(false) {
    fNodes = new Cube::ReconNodeContainerImpl<Cube::ClusterState>;

    // Copy the nodes.  Create new nodes with Cube::ClusterState's
    Cube::ReconNodeContainer::const_iterator in;
    for (in=cluster.GetNodes().begin(); in!=cluster.GetNodes().end();in++){
        Cube::Handle<Cube::ReconNode> node(new Cube::ReconNode);
        Cube::Handle<Cube::ReconObject> object = (*in)->GetObject();
        node->SetObject(object);
        Cube::Handle<Cube::ClusterState> tstate = (*in)->GetState();
        if (tstate){
            Cube::Handle<Cube::ReconState> pstate(
                new Cube::ClusterState(*tstate));
            node->SetState(pstate);
        }
        node->SetQuality((*in)->GetQuality());

        fNodes->push_back(node);
    }


    if (cluster.GetState()){
        Cube::Handle<Cube::ClusterState> state = cluster.GetState();
        fState = new Cube::ClusterState(*state);
    }
    else {
        fState = new Cube::ClusterState;
    }
}

Cube::ReconCluster::~ReconCluster() {}

double Cube::ReconCluster::GetEDeposit() const {
    // I'm being a bit pedantic and casting to the base mix-in class.  This
    // could just as well cast to a Cube::ClusterState.
    const Cube::ClusterState* state
        = dynamic_cast<const Cube::ClusterState*>(fState);
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetEDeposit();
}

double Cube::ReconCluster::GetEDepositVariance() const {
    // I'm being a bit pedantic and casting to the base mix-in class.  This
    // could just as well cast to a Cube::ClusterState.
    const Cube::ClusterState* state
        = dynamic_cast<const Cube::ClusterState*>(fState);
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetEDepositVariance();
}

TLorentzVector Cube::ReconCluster::GetPosition() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ClusterState> state = GetState();
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetPosition();
}

TLorentzVector Cube::ReconCluster::GetPositionVariance() const {
    // This is the preferred way to access a state field.
    Cube::Handle<Cube::ClusterState> state = GetState();
    if (!state) throw std::runtime_error("Cluster state missing");
    return state->GetPositionVariance();
}

bool Cube::ReconCluster::IsXCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.X())) return false;
    return true;
}

bool Cube::ReconCluster::IsYCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Y())) return false;
    return true;
}

bool Cube::ReconCluster::IsZCluster() const {
    TLorentzVector var = GetPositionVariance();
    if (Cube::CorrValues::IsFree(var.Z())) return false;
    return true;
}

int Cube::ReconCluster::GetDimensions() const{
    TLorentzVector var = GetPositionVariance();
    int dim = 0;
    if (IsXCluster()) ++dim;
    if (IsYCluster()) ++dim;
    if (IsZCluster()) ++dim;
    return dim;
}

const Cube::ReconCluster::MomentMatrix& Cube::ReconCluster::GetMoments() const {
    return fMoments;
}

void Cube::ReconCluster::SetMoments(double xx, double yy, double zz,
                                   double xy, double xz, double yz) {
    fTemporariesInitialized = false;
    fMoments(0,0) = xx;
    fMoments(1,1) = yy;
    fMoments(2,2) = zz;
    fMoments(0,1) = xy;
    fMoments(1,0) = xy;
    fMoments(0,2) = xz;
    fMoments(2,0) = xz;
    fMoments(1,2) = yz;
    fMoments(2,1) = yz;
}

void Cube::ReconCluster::SetMoments(const TMatrixT<double>& moments) {
    if (moments.GetNrows() != fMoments.GetNrows()
        || moments.GetNcols() != fMoments.GetNcols()) {
        throw std::runtime_error("Moment matrix sizes do not match.");
    }
    fTemporariesInitialized = false;
    for (int row=0; row<3; ++row) {
        for (int col=0; col<3; ++col) {
            fMoments(row,col) = moments(row,col);
        }
    }
}

void Cube::ReconCluster::UpdateFromHits() {
    fTemporariesInitialized = false;
    // Make sure there is a hit container.
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    if (!hits) return;

    // Make sure the hit container isn't empty.
    Cube::HitSelection::const_iterator begin = hits->begin();
    Cube::HitSelection::const_iterator end = hits->end();
    if (end - begin < 1) return;

    fStatus = Cube::ReconObject::kSuccess;
    fQuality = 1.0;
    fNDOF = std::max(1,int(end - begin - 1));
    Cube::Handle<Cube::ClusterState> state = GetState();

    int dim = state->GetDimensions();
    TVectorT<double> values(dim);
    TVectorT<double> sigmas(dim);

    // Save the index into the state for each of the values.
    int idxE = state->GetEDepositIndex();
    int idxX = state->GetXIndex();
    int idxY = state->GetYIndex();
    int idxZ = state->GetZIndex();
    int idxT = state->GetTIndex();

    // Calculate the weighted mean of the energy deposit, time and position.
    // The stateWeights are chosen as the inverts of the square of the variance
    TVectorT<double> stateValues(dim);
    TVectorT<double> stateNorms(dim);
    for (auto hit = begin; hit != end; ++hit) {
      values(idxE) = (*hit)->GetCharge();
      sigmas(idxE) = 1.; // use unit weight for the charge

      values(idxX) = (*hit)->GetPosition().X();
      sigmas(idxX) = (*hit)->GetUncertainty().X();

      values(idxY) = (*hit)->GetPosition().Y();
      sigmas(idxY) = (*hit)->GetUncertainty().Y();

      values(idxZ) = (*hit)->GetPosition().Z();
      sigmas(idxZ) = (*hit)->GetUncertainty().Z();

      values(idxT) = (*hit)->GetTime();
      sigmas(idxT) = (*hit)->GetTimeUncertainty();

      for (int i = 0; i < dim; ++i) {
        stateValues(i) += values(i) / (sigmas(i) * sigmas(i));
        stateNorms(i) += 1 / (sigmas(i) * sigmas(i));
      }
    }

    for (int i = 0; i < dim; ++i) {
      if (stateNorms(i) > 0)
        stateValues(i) /= stateNorms(i);
    }

    // Estimate the covariance of cluster state values.
    TMatrixTSym<double> stateCov(dim);
    // This counts the weight for each bin of the covariance.
    TMatrixTSym<double> stateWeights(dim);
    // This counts the weight squared for each bin of the covariance.
    TMatrixTSym<double> stateWeights2(dim);

    for (auto hit = begin; hit != end; ++hit) {
      values(idxE) = (*hit)->GetCharge() - stateValues(idxE);
      // The charge deposit is assumed to follow Poisson statistics
      sigmas(idxE) = std::sqrt((*hit)->GetCharge());
      // Treat the charge as known with infinite precision

      values(idxX) = (*hit)->GetPosition().X() - stateValues(idxX);
      sigmas(idxX) = (*hit)->GetUncertainty().X();

      values(idxY) = (*hit)->GetPosition().Y() - stateValues(idxY);
      sigmas(idxY) = (*hit)->GetUncertainty().Y();

      values(idxZ) = (*hit)->GetPosition().Z() - stateValues(idxZ);
      sigmas(idxZ) = (*hit)->GetUncertainty().Z();

      values(idxT) = (*hit)->GetTime() - stateValues(idxT);
      sigmas(idxT) = (*hit)->GetTimeUncertainty();

      for (int row = 0; row < dim; ++row) {
        for (int col = row; col < dim; ++col) {
          double weight = 1.0 / (sigmas(row) * sigmas(col));
          stateCov(row, col) += weight * values(row) * values(col);
          stateWeights(row, col) += weight;
          stateWeights2(row, col) += std::pow(weight, 2);
        }
      }
    }

    // Turn the "stateCov" variable into the sample covariance
    for (int row = 0; row < dim; ++row) {
      for (int col = row; col < dim; ++col) {
        if (row == idxE || col == idxE) {
          stateCov(row, col) = 0.;
        } else {
          stateCov(row, col) *= stateWeights(row, col) /
                                (std::pow(stateWeights(row, col), 2) - stateWeights2(row, col));
        }
        if (col != row)
          stateCov(col, row) = stateCov(row, col);
      }
    }

    // Fix the variance of the deposited energy.  This assumes it's Poisson
    // distributed.
    stateCov(idxE, idxE) = stateValues(idxE);

    // Set the state value and covariance.
    for (int row = 0; row < dim; ++row) {
      state->SetValue(row, stateValues(row));
      for (int col = 0; col < dim; ++col) {
        state->SetCovarianceValue(row, col, stateCov(row, col));
        state->SetCovarianceValue(col, row, stateCov(row, col));
      }
    }

    // Find the moments of the cluster.  This could be done as part of the
    // covariance calculation, but is seperated so it's easier to see what is
    // going on.
    TMatrixTSym<double> moments(3);
    TMatrixTSym<double> chargeSum(3);
    TVector3 center(state->GetValue(idxX),
                    state->GetValue(idxY),
                    state->GetValue(idxZ));

    // Calculate the charge weighted "sum of squared differences".
    for (Cube::HitSelection::const_iterator h = begin;
         h != end;
         ++h) {
        TVector3 diff = (*h)->GetPosition() - center;
        for (int row = 0; row<3; ++row) {
            for (int col = row; col<3; ++col) {
                moments(row,col) += diff(row)*diff(col)*(*h)->GetCharge();
                // Correct for the finite size of the hit.
                if (row==col) {
                    double r = (*h)->GetSize()[row];
                    moments(row,col) += r*r*(*h)->GetCharge();
                }
                chargeSum(row,col) += (*h)->GetCharge();
            }
        }
    }

    // Turn this into the average square difference (i.e. the moments)
    for (int row = 0; row<3; ++row) {
        for (int col = row; col<3; ++col) {
            // No divide by zero please!
            if (chargeSum(row,col) > 1E-6) {
                moments(row,col) = moments(row,col)/chargeSum(row,col);
            }
            else if (row == col) {
                // There were no measurements on this axis so spread the
                // charge over the entire range of the detector
                // (i.e. +-"10 m + epsilon")
                moments(row,col) = 1E+9;
            }
            else {
                // There were no measurements so the aren't any off-axis
                // correlations.
                moments(row,col) = 0;
            }
            // Keep it symmetric.
            moments(col,row) = moments(row,col);
        }
    }

    SetMoments(moments);
}

void Cube::ReconCluster::ls(Option_t *opt) const {
    Cube::ReconObject::ls_base(opt);
    TROOT::IncreaseDirLevel();
    TROOT::IndentLevel();

    if (fState) {
        TROOT::IncreaseDirLevel();
        fState->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    std::ios::fmtflags save = std::cout.flags();
    for (int i = 0; i<3; ++i) {
        TROOT::IndentLevel();
        if (i == 0) std::cout << "Position Cov: ";
        else        std::cout << "              ";
        Cube::Handle<Cube::ClusterState> state = GetState();
        for (int j = 0; j<3; ++j) {
            std::cout << std::setw(12) << state->GetPositionCovariance(i,j);
        }
        std::cout << std::setw(0) << std::endl;
    }
    for (int i = 0; i<3; ++i) {
        TROOT::IndentLevel();
        if (i == 0) std::cout << "Moments:      ";
        else        std::cout << "              ";
        for (int j = 0; j<3; ++j) {
            std::cout << std::setw(12) << fMoments(i,j);
        }
        std::cout << std::setw(0) << std::endl;
    }
    std::cout.flags(save);

    if (fNodes) {
        TROOT::IncreaseDirLevel();
        fNodes->ls(opt);
        TROOT::DecreaseDirLevel();
    }

    TROOT::DecreaseDirLevel();
}

void Cube::ReconCluster::FillTemporaries() const {
    if (fTemporariesInitialized) return;
    fTemporariesInitialized = true;

    TVectorF eigenValues;
    TMatrixF eigenVectors(GetMoments().EigenVectors(eigenValues));

    // Long axis comes first.
    fLongAxis.SetXYZ(eigenVectors(0,0), eigenVectors(1,0), eigenVectors(2,0));
    fLongAxis = std::sqrt(eigenValues(0))*fLongAxis;

    // Major axis comes second.
    fMajorAxis.SetXYZ(eigenVectors(0,1), eigenVectors(1,1), eigenVectors(2,1));
    fMajorAxis = std::sqrt(eigenValues(1))*fMajorAxis;

    // Minor axis comes last.
    fMinorAxis.SetXYZ(eigenVectors(0,2), eigenVectors(1,2), eigenVectors(2,2));
    fMinorAxis = std::sqrt(eigenValues(2))*fMinorAxis;

    const double epsilon = 1E-6;
    if (fLongAxis.X() < -epsilon) {
        fLongAxis = -fLongAxis;
    }
    else if (fLongAxis.X() < epsilon && fLongAxis.Y() < -epsilon) {
        fLongAxis = -fLongAxis;
    }
    else if (fLongAxis.Y() < epsilon && fLongAxis.Z() < 0.0) {
        fLongAxis = -fLongAxis;
    }

    if (fMajorAxis.Y() < -epsilon) {
        fMajorAxis = -fMajorAxis;
    }
    else if (fMajorAxis.Y() < epsilon && fMajorAxis.Z() < -epsilon) {
        fMajorAxis = -fMajorAxis;
    }
    else if (fMajorAxis.Z() < epsilon && fMajorAxis.X() < 0.0) {
        fMajorAxis = -fMajorAxis;
    }

    TVector3 temp = fLongAxis.Cross(fMajorAxis);
    if (temp*fMinorAxis < 0) {
        fMinorAxis = -fMinorAxis;
    }

}

const TVector3& Cube::ReconCluster::GetLongAxis() const {
    FillTemporaries();
    return fLongAxis;
}

const TVector3& Cube::ReconCluster::GetMajorAxis() const {
    FillTemporaries();
    return fMajorAxis;
}

const TVector3& Cube::ReconCluster::GetMinorAxis() const {
    FillTemporaries();
    return fMinorAxis;
}

double Cube::ReconCluster::GetLongExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fLongAxis));
    }
    maxLen /= fLongAxis.Mag();
    return maxLen;

}

double Cube::ReconCluster::GetMajorExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fMajorAxis));
    }
    maxLen /= fMajorAxis.Mag();
    return maxLen;

}

double Cube::ReconCluster::GetMinorExtent() const {
    FillTemporaries();

    double maxLen = 0.0;
    Cube::Handle<Cube::HitSelection> hits = GetHitSelection();
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - GetPosition().Vect();
        maxLen = std::max(maxLen,
                          (*h)->GetSize().Mag() + std::abs(diff*fMinorAxis));
    }
    maxLen /= fMinorAxis.Mag();
    return maxLen;

}
