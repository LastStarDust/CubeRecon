#ifndef G4Hit_hxx_seen
#define G4Hit_hxx_seen

#include "CubeHandle.hxx"

#include <TObject.h>
#include <TLorentzVector.h>

#include <iostream>
#include <vector>

namespace Cube {
    class G4Hit;
}

/// A calibrated hit detector element where the hit information has been
/// reconstructed from the raw hits.  A Hit object represents a hit that has
/// been generated from one or more hits where the position and other
/// properties are calculated from the constituient hits.  A common usage is
/// to represent a hit that has been constructed from two or more hits on a
/// wires that cross at a single X/Y positon.  This is a lighter object that
/// TReconCluster and is primarily used for building up 2D hits into 3D
/// objects.  Complex combinations of hits should be placed into a
/// TReconCluster.
///
/// The Hit class can't be directly instantiated.  It is created
/// using the WritableHit class and is accessed as a Hit class.
class Cube::G4Hit : public TObject {
public:
    G4Hit():
        fSegmentId(-1), fPrimaryId(0), fPrimaryPDG(0),
        fEnergyDeposit(0), fSecondaryDeposit(0),
        fTrackLength(0), fStart(0,0,0,0), fStop(0,0,0,0) {}
    virtual ~G4Hit();

    int GetSegmentId() const {return fSegmentId;}
    void SetSegmentId(int i) {fSegmentId = i;}

    /// The track id of the most important particle associated with this hit
    /// segment.
    int GetPrimaryId() const {return fPrimaryId;}
    void SetPrimaryId(int i) {fPrimaryId = i;}

    /// The main PID
    int GetPDG() const {return fPrimaryPDG;}
    void SetPDG(int i) {fPrimaryPDG = i;}

    /// The total energy deposit in this hit.
    double GetEnergyDeposit() const {return fEnergyDeposit;}
    void SetEnergyDeposit(double v) {fEnergyDeposit = v;}

    /// The "secondary" energy deposit in this hit. Generally, this is used to
    /// help simulate the amount of energy emitted as scintillation light,
    /// i.e. opticalphotons, and is part of the total energy deposit.  The
    /// remaining energy will be deposited as ionization.  In this model (in
    /// argon), the mean number of quanta created will be <N_q> =
    /// (fEnergyDeposit)/(19.5*eV), N_q should be fluctuated around <N_q>,
    /// N_ph = N_q*fSecondaryDeposit/fEnergyDeposit, and N_e = N_q - N_ph.
    /// Thd fSecondaryDeposit value already includes the binomial fluctuation,
    /// so don't fluctuate N_ph or N_e.
    double GetSecondaryDeposit() const {return fSecondaryDeposit;}
    void SetSecondaryDeposit(double v) {fSecondaryDeposit = v;}

    /// The total charged track length in this hit.  This includes the
    /// contribution from all of the secondary particles (e.g. delta-rays)
    /// that are included in this hit.
    double GetTrackLength() const {return fTrackLength;}
    void SetTrackLength(double v) {fTrackLength = v;}

    /// The starting position of the segment.
    const TLorentzVector& GetStart() const {return fStart;}
    void SetStart(const TLorentzVector& v)  {fStart = v;}

    /// The stopping position of the segment.
    const TLorentzVector& GetStop() const {return fStop;}
    void SetStop(const TLorentzVector& v) {fStop = v;}

private:
    /// The segment number.
    Int_t fSegmentId;

    /// The track id of the most important particle associated with this hit
    /// segment.
    Int_t fPrimaryId;

    /// The particle that made the segment.
    Int_t fPrimaryPDG;

    /// The total energy deposit in this hit.
    Float_t fEnergyDeposit;

    /// The "secondary" energy deposit in this hit. Generally, this is used to
    /// help simulate the amount of energy emitted as scintillation light,
    /// i.e. opticalphotons, and is part of the total energy deposit.  The
    /// remaining energy will be deposited as ionization.  In this model (in
    /// argon), the mean number of quanta created will be <N_q> =
    /// (fEnergyDeposit)/(19.5*eV), N_q should be fluctuated around <N_q>,
    /// N_ph = N_q*fSecondaryDeposit/fEnergyDeposit, and N_e = N_q - N_ph.
    /// Thd fSecondaryDeposit value already includes the binomial fluctuation,
    /// so don't fluctuate N_ph or N_e.
    Float_t fSecondaryDeposit;

    /// The total charged track length in this hit.  This includes the
    /// contribution from all of the secondary particles (e.g. delta-rays)
    /// that are included in this hit.
    Float_t fTrackLength;

    /// The starting position of the segment.
    TLorentzVector fStart;

    /// The stopping position of the segment.
    TLorentzVector fStop;

    ClassDef(G4Hit,1);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
