#ifndef CubeReconState_hxx_seen
#define CubeReconState_hxx_seen
////////////////////////////////////////////////////////////////
// A class to hold the state information for the reconstruction objects.  This
// file defines several classes and needs to be split into multiple files
// after the design is complete.
#include <vector>
#include <string>

#include <TObject.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "CubeCorrValues.hxx"
#include "CubeHandle.hxx"

namespace Cube {
    class ReconState;
    class ClusterState;
    class ShowerState;
    class TrackState;
    class PIDState;
    class VertexState;
}

/// The ReconState class represents the value of parameters and covariances
/// associated with reconstruction objects.  It does not contain general
/// status information such as the goodness, degree's of freedom, or objects
/// which contributed to the reconstruction.  The values saved in the particle
/// states are described in \ref recObjGeometric and \ref recObjPhysics.
///
/// This is the polymorphic base class that all other state classes are
/// derived from, and it provides minimal operations.  The main purpose of
/// this class is to allow polymorphic vectors of states.  However, it
/// provides some minimal access to the contained data.
class Cube::ReconState: public TObject {
public:
    ReconState();
    ReconState(const Cube::ReconState& state);
    virtual ~ReconState();

    /// Return a string with all of the state field names.  This name is used
    /// to build a type has for the state.
    std::string GetStateFields() const;

    /// Return the number of dimensions in the state (the size of the state).
    int GetDimensions() const;

    /// Not the preferred interface, but get the value by index.  The best way
    /// to get a value is through the Get method (e.g. GetEDeposit());
    double GetValue(int i) const;

    /// Not the preferred interface, but set the value by index.  The best way
    /// to get a value is through the Set method provided by the mix-in
    /// classes (e.g. SetEDeposit()).
    void SetValue(int i, double value);

    /// Not the preferred interface: get the covariance by index.
    double GetCovarianceValue(int i, int j) const;

    /// Not the preferred interface: set the covariance by index.
    void SetCovarianceValue(int i, int j, double value);

    /// The projection operator to get the full state.
    static Cube::CorrValues ProjectState(
        const Cube::Handle<ReconState>& state);

    /// Set a parameter to be free (unconstrained).
    void SetFree(int i);

    /// Check if a parameter is free.
    bool IsFree(int i) const;

    /// Check if a variance corresponds to a free parameter.
    bool IsFree(double v) const;

    /// Set a parameter to be fixed.
    void SetFixed(int i);

    /// Check if a parameter is fixed.
    bool IsFixed(int i) const;

    /// Check if a variance corresponds to a fixed parameter.
    bool IsFixed(double v) const;

    /// Validate the covariance.
    void Validate();

    /// Print the object information.
    virtual void ls(Option_t *opt = "") const;
    virtual void Print(Option_t *opt = "") const {ls(opt);}

protected:
    friend class ClusterState;
    friend class ShowerState;
    friend class TrackState;
    friend class PIDState;
    friend class VertexState;

    /// A final initialization routine that is called in the constructor of
    /// the instantiated class.  This builds the actual state vector.
    void Init();

    /// The vector of correlated values (a vector and a covariance) that holds
    /// the state information.
    Cube::CorrValues fValues;

    /// A vector of parameter names.  This identifies the fields in the state.
    std::vector<std::string> fFieldNames;

    ClassDef(ReconState,1);
};

/// A macro that adds the energy deposit property to a state.
/// The position state is defined as a scaler and the
/// associated uncertainty.
///
/// \bug The units of the energy deposit depend on the context.  After
/// particle ID, the deposit can be in uncalibrated or calibrated energy.  In
/// the TPC the deposit might be in terms of fempto-coulombs, but for a PMT
/// the deposit might be based on photo-electrons.
#define ENERGY_DEPOSIT_STATE_DECLARATION                                \
    public:                                                             \
    /**  Get the index of the EDeposit field in the CorrValues vector.*/ \
    int GetEDepositIndex() const {return fEDepositIndex;}               \
    /** Get the value of the EDeposit. */                               \
    double GetEDeposit() const {return GetValue(fEDepositIndex);}       \
    /** Set the value of the EDeposit. */                               \
    void SetEDeposit(double enr) {SetValue(fEDepositIndex,enr);}        \
    /** Get the variance of the EDeposit.*/                             \
    double GetEDepositVariance(void) const                              \
    {return GetCovarianceValue(fEDepositIndex,fEDepositIndex);}         \
    /** Set the variance of the EDeposit.  */                           \
    void SetEDepositVariance(double var)                                \
    {SetCovarianceValue(fEDepositIndex,fEDepositIndex,var);}

#define ENERGY_DEPOSIT_STATE_PRIVATE                                    \
    /** The index of the EDeposit field in the CorrValues vector.*/    \
    private: unsigned char fEDepositIndex

/// This should be included in the class constructor.
#define ENERGY_DEPOSIT_STATE_DEFINITION                                 \
    fEDepositIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("EDeposit")

/// A macro that adds the position property to a state.  The position state is
/// defined as a four vector (X,Y,Z,T) and the associated uncertainties.
#define POSITION_STATE_DECLARATION                                      \
    public:                                                             \
    /** Get the size of the position field */                           \
    int GetPositionSize() const {return 4;}                             \
    /** Get the index of the position field.*/                          \
    int GetPositionIndex() const {return fPositionIndex;}               \
    /** Get the index of the X,Y,Z,T fields in the CorrValues vector.*/ \
    int GetXIndex() const {return fPositionIndex;}                      \
    int GetYIndex() const {return fPositionIndex+1;}                    \
    int GetZIndex() const {return fPositionIndex+2;}                    \
    int GetTIndex() const {return fPositionIndex+3;}                    \
    /** Get the value of the Position. */                               \
    TLorentzVector GetPosition() const {                                \
        return TLorentzVector(GetValue(fPositionIndex+0),               \
                              GetValue(fPositionIndex+1),               \
                              GetValue(fPositionIndex+2),               \
                              GetValue(fPositionIndex+3));              \
    }                                                                   \
    /** Set the value of the Position.*/                                \
    void SetPosition(double x, double y, double z, double t) {          \
        SetValue(fPositionIndex+0,x);                                   \
        SetValue(fPositionIndex+1,y);                                   \
        SetValue(fPositionIndex+2,z);                                   \
        SetValue(fPositionIndex+3,t);                                   \
    }                                                                   \
    /** Set the value of the Position.*/                                \
    void SetPosition(const TLorentzVector& pos) {                       \
        SetPosition(pos.X(),pos.Y(),pos.Z(),pos.T());                   \
    }                                                                   \
    /** Set the covariance of the position variables.  The indices */   \
    /* run from 0 to 3 (0: X, 1: Y, 2: Z, 3: T).*/                      \
    void SetPositionCovariance(int i, int j, double v) {                \
        SetCovarianceValue(fPositionIndex+i,fPositionIndex+j,v);        \
    }                                                                   \
    /** Get the covariance of the position variables.  The indices */   \
    /* run from 0 to 3 (0: X, 1: Y, 2: Z, 3: T).*/                      \
    double GetPositionCovariance(int i, int j) const {                  \
        return GetCovarianceValue(fPositionIndex+i,fPositionIndex+j);   \
    }                                                                   \
    /** Get the variance of the Position.*/                             \
    TLorentzVector GetPositionVariance(void) const {                    \
        return TLorentzVector(                                          \
            GetPositionCovariance(0,0),                                 \
            GetPositionCovariance(1,1),                                 \
            GetPositionCovariance(2,2),                                 \
            GetPositionCovariance(3,3));                                \
    }                                                                   \
    /** Set the variance of the Position.*/                             \
    void SetPositionVariance(double x, double y, double z, double t) {  \
        SetPositionCovariance(0,0,x);                                   \
        SetPositionCovariance(1,1,y);                                   \
        SetPositionCovariance(2,2,z);                                   \
        SetPositionCovariance(3,3,t);                                   \
    }

#define POSITION_STATE_PRIVATE                                          \
    private: unsigned char fPositionIndex

/// This should be included in the class constructor.
#define POSITION_STATE_DEFINITION                                       \
    fPositionIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("X");                                         \
    fFieldNames.push_back("Y");                                         \
    fFieldNames.push_back("Z");                                         \
    fFieldNames.push_back("T")

/// A macro that adds the direction property to a state.  The direction state
/// is defined as a three vector (dX,dY,dZ) and the associated uncertainties.
/// The correlations between the components are enforced in the covariance
/// matrix.
#define DIRECTION_STATE_DECLARATION                                     \
    public:                                                             \
    /** Get the size of the direction field */                          \
    int GetDirectionSize() const {return 3;}                            \
    /** Get the index of the direction field.*/                         \
    int GetDirectionIndex() const {return fDirectionIndex;}             \
    /** Get the value of the direction. */                              \
    TVector3 GetDirection() const {                                     \
        return TVector3(GetValue(fDirectionIndex+0),                    \
                        GetValue(fDirectionIndex+1),                    \
                        GetValue(fDirectionIndex+2));                   \
    }                                                                   \
    /** Set the value of the direction.*/                               \
    void SetDirection(double x, double y, double z) {                   \
        SetValue(fDirectionIndex+0,x);                                  \
        SetValue(fDirectionIndex+1,y);                                  \
        SetValue(fDirectionIndex+2,z);                                  \
    }                                                                   \
    /** Set the value of the direction.*/                               \
    void SetDirection(const TVector3& dir) {                            \
        SetDirection(dir.X(),dir.Y(),dir.Z());                          \
    }                                                                   \
    /** Set the covariance of the direction variables.  The indices */  \
    /* run from 0 to 2 (0: dX, 1: dY, 2: dZ).*/                         \
    void SetDirectionCovariance(int i, int j, double v) {               \
        SetCovarianceValue(fDirectionIndex+i,fDirectionIndex+j,v);      \
    }                                                                   \
    /** Get the covariance of the direction variables.  The indices */  \
    /* run from 0 to 2 (0: dX, 1: dY, 2: dZ).*/                         \
    double GetDirectionCovariance(int i, int j) const {                 \
        return GetCovarianceValue(fDirectionIndex+i,fDirectionIndex+j); \
    }                                                                   \
    /** Get the variance of the Direction.*/                            \
    TVector3 GetDirectionVariance(void) const {                         \
        return TVector3(                                                \
            GetDirectionCovariance(0,0),                                \
            GetDirectionCovariance(1,1),                                \
            GetDirectionCovariance(2,2));                               \
    }                                                                   \
    /** Set the variance of the Direction.*/                            \
    void SetDirectionVariance(double x, double y, double z) {           \
        SetDirectionCovariance(0,0,x);                                  \
        SetDirectionCovariance(1,1,y);                                  \
        SetDirectionCovariance(2,2,z);                                  \
    }

#define DIRECTION_STATE_PRIVATE                 \
    private: unsigned char fDirectionIndex

/// This should be included in the class constructor.
#define DIRECTION_STATE_DEFINITION                                      \
    fDirectionIndex=fFieldNames.size();                                 \
    fFieldNames.push_back("dX");                                        \
    fFieldNames.push_back("dY");                                        \
    fFieldNames.push_back("dZ")

/// A macro that adds the width of a TReconShower object.  The cone value
/// depends on the type of the shower fit.  For a EM fit, the cone value will
/// represent the local width of the shower.
#define CONE_STATE_DECLARATION                                          \
    public:                                                             \
    /** Get the size of the cone field */                               \
    int GetConeSize() const {return 1;}                                 \
    /** Get the index of the cone field in the CorrValues vector.*/    \
    int GetConeIndex() const {return fConeIndex;}                       \
    /** Get the value of the cone. */                                   \
    double GetCone() const {return GetValue(fConeIndex);}               \
    /** Set the value of the cone. */                                   \
    void SetCone(double enr) {SetValue(fConeIndex,enr);}                \
    /** Get the variance of the cone.*/                                 \
    double GetConeVariance(void) const                                  \
    {return GetCovarianceValue(fConeIndex,fConeIndex);}                 \
    /** Set the variance of the cone.  */                               \
    void SetConeVariance(double var)                                    \
    {SetCovarianceValue(fConeIndex,fConeIndex,var);}

#define CONE_STATE_PRIVATE                      \
    private: unsigned char fConeIndex

/// This should be included in the class constructor.
#define CONE_STATE_DEFINITION                                       \
    fConeIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("Cone")

/// A macro that adds the mass to a state.  The mass state also holds the
/// associated uncertainties.
#define MASS_STATE_DECLARATION                                          \
    public:                                                             \
    /** Get the size of the mass field */                               \
    int GetMassSize() const {return 1;}                                 \
    /** Get the index of the mass field in the CorrValues vector.*/     \
    int GetMassIndex() const {return fMassIndex;}                       \
    /** Get the value of the mass. */                                   \
    double GetMass() const {return GetValue(fMassIndex);}               \
    /** Set the value of the mass. */                                   \
    void SetMass(double mass) {SetValue(fMassIndex,mass);}              \
    /** Get the variance of the mass. */                                \
    double GetMassVariance(void) const                                  \
    {return GetCovarianceValue(fMassIndex,fMassIndex);}                 \
    /** Set the variance of the mass. */                                \
    void SetMassVariance(double mass)                                   \
    {SetCovarianceValue(fMassIndex,fMassIndex,mass);}

#define MASS_STATE_PRIVATE                      \
    private: unsigned char fMassIndex

/// This should be included in the class constructor.
#define MASS_STATE_DEFINITION                                       \
    fMassIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("Mass")

/// A macro that adds the curvature property to a state.  The curvature state
/// is defined as a three vector (cX,cY,cZ) and the associated uncertainties.
/// The correlations between the components are enforced in the covariance
/// matrix.
#define CURVATURE_STATE_DECLARATION                                     \
    public:                                                             \
    /** Get the size of the direction field */                          \
    int GetCurvatureSize() const {return 3;}                            \
    /** Get the index of the direction field.*/                         \
    int GetCurvatureIndex() const {return fCurvatureIndex;}             \
    /** Get the value of the direction. */                              \
    TVector3 GetCurvature() const {                                     \
        return TVector3(GetValue(fCurvatureIndex+0),                    \
                        GetValue(fCurvatureIndex+1),                    \
                        GetValue(fCurvatureIndex+2));                   \
    }                                                                   \
    /** Set the value of the direction.*/                               \
    void SetCurvature(double x, double y, double z) {                   \
        SetValue(fCurvatureIndex+0,x);                                  \
        SetValue(fCurvatureIndex+1,y);                                  \
        SetValue(fCurvatureIndex+2,z);                                  \
    }                                                                   \
    /** Set the value of the direction.*/                               \
    void SetCurvature(const TVector3& curv) {                           \
        SetCurvature(curv.X(),curv.Y(),curv.Z());                       \
    }                                                                   \
    /** Set the covariance of the direction variables.  The indices */  \
    /* run from 0 to 2 (0: dX, 1: dY, 2: dZ).*/                         \
    void SetCurvatureCovariance(int i, int j, double v) {               \
        SetCovarianceValue(fCurvatureIndex+i,fCurvatureIndex+j,v);      \
    }                                                                   \
    /** Get the covariance of the direction variables.  The indices */  \
    /* run from 0 to 2 (0: dX, 1: dY, 2: dZ).*/                         \
    double GetCurvatureCovariance(int i, int j) const {                 \
        return GetCovarianceValue(fCurvatureIndex+i,fCurvatureIndex+j); \
    }                                                                   \
    /** Get the variance of the Curvature.*/                            \
    TVector3 GetCurvatureVariance(void) const {                         \
        return TVector3(                                                \
            GetCurvatureCovariance(0,0),                                \
            GetCurvatureCovariance(1,1),                                \
            GetCurvatureCovariance(2,2));                               \
    }                                                                   \
    /** Set the variance of the Curvature.*/                            \
    void SetCurvatureVariance(double x, double y, double z) {           \
        SetCurvatureCovariance(0,0,x);                                  \
        SetCurvatureCovariance(1,1,y);                                  \
        SetCurvatureCovariance(2,2,z);                                  \
    }

#define CURVATURE_STATE_PRIVATE                 \
    private: unsigned char fCurvatureIndex

/// This should be included in the class constructor.
#define CURVATURE_STATE_DEFINITION                                      \
    fCurvatureIndex=fFieldNames.size();                                 \
    fFieldNames.push_back("curvX");                                     \
    fFieldNames.push_back("curvY");                                     \
    fFieldNames.push_back("curvZ")

/// A macro that adds the width of a curvilinear energy deposit property to a
/// state.  The width is the extent of an energy deposition perpendicular to
/// local direction and the associated uncertainties.
#define WIDTH_STATE_DECLARATION                                         \
    public:                                                             \
    /** Get the size of the Position field */                           \
    int GetWidthSize() const {return 1;}                                \
    /** Get the index of the width field in the TCorrValuest vector.*/  \
    int GetWidthIndex() const {return fWidthIndex;}                     \
    /** Get the value of the width. */                                  \
    double GetWidth() const {return GetValue(fWidthIndex);}             \
    /** Set the value of the width. */                                  \
    void SetWidth(double width) {SetValue(fWidthIndex,width);}          \
    /** Get the variance of the width. */                               \
    double GetWidthVariance(void) const                                 \
    {return GetCovarianceValue(fWidthIndex,fWidthIndex);}               \
    /** Set the variance of the width. */                               \
    void SetWidthVariance(double width)                                 \
    {SetCovarianceValue(fWidthIndex,fWidthIndex,width);}

#define WIDTH_STATE_PRIVATE                     \
    private: unsigned char fWidthIndex

/// This should be included in the class constructor.
#define WIDTH_STATE_DEFINITION                                       \
    fWidthIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("Width")

/// A macro that adds a property for the magnitude of the momentum to a state.
#define MOMENTUM_STATE_DECLARATION                                      \
    public:                                                             \
    /** Get the size of the momentum field */                           \
    int GetMomentumSize() const {return 1;}                             \
    /** Get the index of the momentum field in the TCorrValues vector.*/ \
    int GetMomentumIndex() const {return fMomentumIndex;}               \
    /** Get the value of the momentum. */                               \
    double GetMomentum() const {return GetValue(fMomentumIndex);}       \
    /** Set the value of the momentum. */                               \
    void SetMomentum(double momentum) {SetValue(fMomentumIndex,momentum);} \
    /** Get the variance of the momentum. */                            \
    double GetMomentumVariance(void) const                              \
    {return GetCovarianceValue(fMomentumIndex,fMomentumIndex);}         \
    /** Set the variance of the momentum. */                            \
    void SetMomentumVariance(double momentum)                           \
    {SetCovarianceValue(fMomentumIndex,fMomentumIndex,momentum);}

#define MOMENTUM_STATE_PRIVATE                                          \
    private: unsigned char fMomentumIndex

/// This should be included in the class constructor.
#define MOMENTUM_STATE_DEFINITION                                       \
    fMomentumIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("Momentum")

/// A macro that adds a property for the magnitude of the particle charge to a
/// state.
#define CHARGE_STATE_DECLARATION                                        \
    public:                                                             \
    /** Get the size of the charge field */                             \
    int GetChargeSize() const {return 4;}                               \
    /** Get the index of the charge field in the TCorrValues vector.*/  \
    int GetChargeIndex() const {return fChargeIndex;}                   \
    /** Get the value of the charge. */                                 \
    double GetCharge() const {return GetValue(fChargeIndex);}           \
    /** Set the value of the charge. */                                 \
    void SetCharge(double charge) {SetValue(fChargeIndex,charge);}      \
    /** Get the variance of the charge. */                              \
    double GetChargeVariance(void) const                                \
    {return GetCovarianceValue(fChargeIndex,fChargeIndex);}             \
    /** Set the variance of the charge. */                              \
    void SetChargeVariance(double charge)                               \
    {SetCovarianceValue(fChargeIndex,fChargeIndex,charge);}

#define CHARGE_STATE_PRIVATE                                           \
    private: unsigned char fChargeIndex

/// This should be included in the class constructor.
#define CHARGE_STATE_DEFINITION                                       \
    fChargeIndex=fFieldNames.size();                                  \
    fFieldNames.push_back("Charge")

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
