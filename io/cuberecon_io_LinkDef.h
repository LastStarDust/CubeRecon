#ifdef __CINT__
#include "CubeHandle.hxx"
#include "CubeEvent.hxx"
#include "CubeAlgorithm.hxx"
#include "CubeAlgorithmResult.hxx"
#include "CubeHit.hxx"
#include "CubeHitSelection.hxx"
#include "CubeCorrValues.hxx"
#include "CubeReconState.hxx"
#include "CubeVertexState.hxx"
#include "CubeClusterState.hxx"
#include "CubeShowerState.hxx"
#include "CubeTrackState.hxx"
#include "CubeReconObject.hxx"
#include "CubeReconNode.hxx"
#include "CubeReconVertex.hxx"
#include "CubeReconCluster.hxx"
#include "CubeReconShower.hxx"
#include "CubeReconTrack.hxx"
#include "CubeG4Hit.hxx"
#include "CubeG4Trajectory.hxx"

#pragma link C++ class Cube::HandleBase+;
#pragma link C++ class Cube::HandleBaseDeletable+;
#pragma link C++ class Cube::HandleBaseUndeletable+;
#pragma link C++ class Cube::VHandle+;

#pragma link C++ class Cube::Event;
#pragma link C++ class Cube::Handle<Cube::Event>;

#pragma link C++ class Cube::Algorithm;
#pragma link C++ class Cube::Handle<Cube::Algorithm>;

#pragma link C++ class Cube::AlgorithmResult+;
#pragma link C++ class Cube::Handle<Cube::AlgorithmResult>+;

#pragma link C++ class Cube::Hit+;
#pragma link C++ class Cube::Handle<Cube::Hit>+;

#pragma link C++ class Cube::WritableHit+;
#pragma link C++ class Cube::Handle<Cube::WritableHit>+;

#pragma link C++ class Cube::HitSelection+;
#pragma link C++ class Cube::Handle<Cube::HitSelection>+;

#pragma link C++ class Cube::CorrValues+;

#pragma link C++ class Cube::ReconState+;
#pragma link C++ class Cube::Handle<Cube::ReconState>+;

#pragma link C++ class Cube::VertexState+;
#pragma link C++ class Cube::Handle<Cube::VertexState>+;
#pragma link C++ class Cube::ReconNodeContainerImpl<Cube::VertexState>+;

#pragma link C++ class Cube::ClusterState+;
#pragma link C++ class Cube::Handle<Cube::ClusterState>+;
#pragma link C++ class Cube::ReconNodeContainerImpl<Cube::ClusterState>+;

#pragma link C++ class Cube::ShowerState+;
#pragma link C++ class Cube::Handle<Cube::ShowerState>+;
#pragma link C++ class Cube::ReconNodeContainerImpl<Cube::ShowerState>+;

#pragma link C++ class Cube::TrackState+;
#pragma link C++ class Cube::Handle<Cube::TrackState>+;
#pragma link C++ class Cube::ReconNodeContainerImpl<Cube::TrackState>+;

#pragma link C++ class Cube::ReconNode+;
#pragma link C++ class Cube::Handle<Cube::ReconNode>+;
#pragma link C++ class Cube::ReconNodeContainer+;

#pragma link C++ class Cube::ReconObject+;
#pragma link C++ class Cube::Handle<Cube::ReconObject>+;

#pragma link C++ class Cube::ReconObjectContainer+;
#pragma link C++ class Cube::Handle<Cube::ReconObjectContainer>+;

#pragma link C++ class Cube::ReconVertex+;
#pragma link C++ class Cube::Handle<Cube::ReconVertex>+;

#pragma link C++ class Cube::ReconCluster+;
#pragma link C++ class Cube::Handle<Cube::ReconCluster>+;

#pragma link C++ class Cube::ReconShower+;
#pragma link C++ class Cube::Handle<Cube::ReconShower>+;

#pragma link C++ class Cube::ReconTrack+;
#pragma link C++ class Cube::Handle<Cube::ReconTrack>+;

#pragma link C++ class Cube::G4Hit+;
#pragma link C++ class Cube::Handle<Cube::G4Hit>+;

#pragma link C++ class Cube::G4Trajectory+;
#pragma link C++ class Cube::Handle<Cube::G4Trajectory>+;

#endif
