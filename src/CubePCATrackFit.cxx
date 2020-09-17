#include "CubePCATrackFit.hxx"

#include <CubeReconNode.hxx>
#include <CubeReconCluster.hxx>
#include <CubeLog.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>

#include <memory>

namespace {

    // Take a position and turn it into a principal component value.
    double FindPrincipal(TPrincipal& pca, const TVector3& position) {
        double X[3] = {position.X(), position.Y(), position.Z()};
        double P[3] = {0,0,0};
        pca.X2P(X,P);
        return P[0];
    }

    // Take a principal component value and turn it into a position.
    TVector3
    FindPosition(TPrincipal& pca, double principal) {
        double X[3] = {0,0,0};
        double P[3] = {principal,0,0};
        pca.P2X(P,X,3);
        return TVector3(X[0],X[1],X[2]);
    }
}

Cube::PCATrackFit::PCATrackFit() {}
Cube::PCATrackFit::~PCATrackFit() {}

Cube::Handle<Cube::ReconTrack>
Cube::PCATrackFit::Apply(Cube::Handle<Cube::ReconTrack>& input) {

    Cube::ReconNodeContainer& nodes = input->GetNodes();
    if (nodes.size() < 2) {
        CUBE_ERROR << "Not enough nodes to fit." << std::endl;
        return Cube::Handle<Cube::ReconTrack>();
    }

    CUBE_LOG(2) << "Start PCA fit with " << nodes.size() << " nodes"
                << std::endl;

    /////////////////////////////////////////////////////////////////////
    /// Fill the PCA using the node object positions.  This also makes a very
    /// crude estimate of the position covariance.
    /////////////////////////////////////////////////////////////////////
    std::unique_ptr<TPrincipal> pca(new TPrincipal(3,""));
    TMatrixD posCov(3,3);
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Get the object from the node.  It had better be a cluster.
        Cube::Handle<Cube::ReconCluster> cluster = (*n)->GetObject();
        if (!cluster) {
            CUBE_ERROR << "Missing cluster in track" << std::endl;
            return Cube::Handle<Cube::ReconTrack>();
        }
        double row1[3] = {cluster->GetPosition().X(),
                          cluster->GetPosition().Y(),
                          cluster->GetPosition().Z()};
        for (double q = cluster->GetEDeposit(); q > 0; q -= 1) {
            pca->AddRow(row1);
        }

        Cube::ReconCluster::MomentMatrix moments = cluster->GetMoments();
        moments.InvertFast();
        posCov += moments;
    }
    pca->MakePrincipals();
    posCov.InvertFast();

    /////////////////////////////////////////////////////////////////////
    /// Use the PCA to find the direction of the track.  This also makes a
    /// crude estimate of the direction covariance.
    /////////////////////////////////////////////////////////////////////

    /// Find the front and last positions of the track by projecting using the
    /// PCA.
    Cube::Handle<Cube::ReconCluster> frontCluster = nodes.front()->GetObject();
    double frontPrincipal
        = FindPrincipal(*pca, frontCluster->GetPosition().Vect());
    TVector3 frontPosition = FindPosition(*pca,frontPrincipal);

    Cube::Handle<Cube::ReconCluster> backCluster = nodes.back()->GetObject();
    double backPrincipal
        = FindPrincipal(*pca, backCluster->GetPosition().Vect());
    TVector3 backPosition = FindPosition(*pca,backPrincipal);

    // The track direction is just the direction between the front and back
    // ends of the track.
    TVector3 dir = (backPosition-frontPosition).Unit();
    TMatrixD dirCov(3,3);
    double length2 = (backPosition-frontPosition).Mag2();
    Cube::Handle<Cube::ClusterState> frontState = frontCluster->GetState();
    Cube::Handle<Cube::ClusterState> backState = backCluster->GetState();
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            double v = frontCluster->GetMoments()(i,j)
                + backCluster->GetMoments()(i,j);
            v = v/length2;
            dirCov(i,j) = v;
        }
    }

    /////////////////////////////////////////////////////////////////////
    /// Use the PCA and direction to fill the state at each node.
    /////////////////////////////////////////////////////////////////////
    double logLikelihood = 0.0;
    double energyDeposit = 0.0;
    double energyVariance = 0.0;
    for (Cube::ReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Get the state to be filled from the node.  It had better be a
        // trackState!.
        Cube::Handle<Cube::TrackState> trackState = (*n)->GetState();

        // Get the object from the node.  It had better be a cluster.
        Cube::Handle<Cube::ReconCluster> cluster = (*n)->GetObject();
        Cube::Handle<Cube::ClusterState> clusterState = cluster->GetState();

        // Find the position for this node.
        double nodePrincipal
            = FindPrincipal(*pca, clusterState->GetPosition().Vect());
        TVector3 nodePosition = FindPosition(*pca,nodePrincipal);

        // Sum the track energy deposit and variance so that we can fill the
        // track state later.
        energyDeposit += clusterState->GetEDeposit();
        energyVariance += clusterState->GetEDepositVariance();

        // Set the track state using the estimated node position and direction
        trackState->SetEDeposit(clusterState->GetEDeposit());
        trackState->SetPosition(nodePosition.X(),
                                nodePosition.Y(),
                                nodePosition.Z(),
                                clusterState->GetPosition().T());
        trackState->SetDirection(dir.X(),dir.Y(),dir.Z());

        // Set the track state covariances.
        trackState->SetEDepositVariance(clusterState->GetEDepositVariance());
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = posCov(i,j);
                trackState->SetPositionCovariance(i,j,v);
                v = dirCov(i,j);
                trackState->SetDirectionCovariance(i,j,v);
            }
        }

        // Calculate the goodness.  This depends on the idea that the cluster
        // covariance is always diagonal.
        TVector3 nodeDiff = nodePosition - cluster->GetPosition().Vect();
        for (int i=0; i<3; ++i) {
            logLikelihood +=
                nodeDiff[i]*nodeDiff[i]/cluster->GetPositionVariance()[i];
        }
    }

    // Fill the overall track state at the front
    Cube::Handle<Cube::TrackState> trackState = input->GetState();
    Cube::Handle<Cube::TrackState> nodeState = nodes.front()->GetState();
    *trackState = *nodeState;
    trackState->SetEDeposit(energyDeposit);
    trackState->SetEDepositVariance(energyVariance);

    // Fill the overall track state at the back
    trackState = input->GetBack();
    nodeState = nodes.back()->GetState();
    *trackState = *nodeState;
    trackState->SetEDeposit(energyDeposit);
    trackState->SetEDepositVariance(energyVariance);

    // Setup the track information and status fields.
    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(Cube::ReconObject::kSuccess);
    input->SetStatus(Cube::ReconObject::kRan);
    input->SetAlgorithmName("PCATrackFit");
    input->SetQuality(logLikelihood);
    input->SetNDOF(trackDOF);

    return input;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
