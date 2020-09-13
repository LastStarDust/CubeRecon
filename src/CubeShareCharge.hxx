#ifndef CubeShareCharge_hxx_seen
#define CubeShareCharge_hxx_seen
/////////////////////////////////////////////////////////////////////
//
// Distribute the charge on each fiber to the cubes.  Takes a HitSelection
// which contains composit hits constructed.  The hits must contain the simple
// Hit objects that were used to construct the hit
//
// Jargon Alert:
//
//   Deposit == The energy deposited in a cube (or the deposited energy going
//       into a fiber).  This is the energy (measured in "pe") that was
//       deposited into a cube by a particle (even if it's leaking through the
//       walls).  The deposit does not include the effect of attenuation.  The
//       deposit going into a fiber is 1/3 the total deposit in a cube.
//
//   Measurement == The charge (measured in "pe") seen by the sensor.  This
//       can either be the expected charge, or the actual measurement made by
//       the sensor.  It should be clear from the context (expected vs
//       observed).
//
/////////////////////////////////////////////////////////////////////
#include <CubeHitSelection.hxx>
#include <CubeHit.hxx>

#include <memory>
#include <vector>
#include <map>

namespace Cube {
    class ShareCharge;
};

class Cube::ShareCharge {
private:
    struct AugmentedCube;
    struct AugmentedDeposit;
    struct AugmentedFiber;

    // This is leftover from the original implementation using globals, and
    // won't work!
    std::vector<std::shared_ptr<AugmentedCube>> fAugmentedCubes;
    std::vector<std::shared_ptr<AugmentedDeposit>> fAugmentedDeposits;
    std::map<Cube::Handle<Cube::Hit>,
             std::shared_ptr<AugmentedFiber>> fAugmentedFibers;

    // An object to keep track of the energy deposited in each cube.  The
    // deposited energy is then split equally between the fibers that are
    // reading out this cube.  The link between the cube and fiber is tracked
    // using the AugmentedDeposit class.
    struct AugmentedCube {
        // The associated Hit for this cube.  The charge for this hit will be
        // adjusted based on the final deposit.
        Cube::Handle<Cube::Hit> Hit;

        // The index of the associated deposits in the vector of
        // AugmentedDeposit objects.
        std::vector<std::weak_ptr<AugmentedDeposit> > Deposits;

        // The derivative of changing *one* deposit based on the other
        // deposits in the cube.  This is actually just a function, but is
        // here since it is conceptually associated with the cube.  If one of
        // the fibers is missing, then the associated charge should be set to
        // zero (or negative).  The formula is generated using maxima, and is
        // the derivative of
        //
        // Qavg : (q1+q2+q3)/3
        // X2: [(q1-Qavg)^2 + (q2-Qavg)^2 + (q3+Qavg)^2]/(Qavg + 1)
        // f90(diff(X2,q1));
        //
        // The code is then tweaked by hand to make sure that it fits C++.  The
        // function artificially imposes the constraint that all of the charges
        // are greater or equal to zero.
        double CubeDepositDerivative(double q1, double q2, double q3) const;

        // Get the total deposit in this cube.  The deposits have units of
        // "pe", but correspond to the *energy* deposited in the cube.  This
        // is the sum of the deposited energy associated with each fiber.
        double GetDeposit() const;

        // Set the total deposit in this cube.  This changes the deposit for
        // each fiber contributing to this cube, and imposes the constraint
        // that each fiber measures the same energy in the cube.  (i.e. if
        // three fibers contribute to this cube, the deposit associated with
        // each fiber will be one third of the total cube deposit).
        void SetDeposit(double dep);
    };

    // An object to keep track of the energy deposits that go to a particular
    // fiber.  Once this is created, the calculation will not be accessing the
    // fiber because the only information that is used is the measured charge
    // in the fiber.  This object keeps track of which deposits contributed to
    // the fiber.
    struct AugmentedFiber {
        // The THit for the fiber.
        Cube::Handle<Cube::Hit> Hit;
        // The measured number of photo electrons in the hit (this is a
        // copy so we don't need to keep accessing the hit.
        double Measurement;
        // A vector of the indexes of the objects in the AugmentedDeposit
        // vector that are connected to cubes on this fiber.
        std::vector<std::weak_ptr<AugmentedDeposit> > Deposits;
        // Check to see if the fiber is shared by multiple cubes.  If the
        // fiber isn't shared, then there is only one deposit associated with
        // this fiber, and it's value should be fixed to this measurement.
        bool IsSharedFiber() const {return Deposits.size()>1;}
        // Get the measurement for this fiber.
        double GetMeasurement() const {return Measurement;}
        // Get the current estimate of the expected measurements from
        // different cubes into this fiber.  When the calculation has
        // converged the sum of the expected measurements will be equal to the
        // measurement.
        double GetExpectedMeasurement() const;
    };

    // An object to keep track of the contribution of the deposited energy
    // from one cube to one particular fiber. Since a fiber may be measuring
    // several cube deposits, it may contributions from several deposite.
    // This is used to track the linkage between the 2D and 3D hits.  The
    // basic information is saved as the energy deposition being contributed
    // by the cube, and the attenuation correction needed to turn this into
    // the contribution to the measurement at the fiber.
    struct AugmentedDeposit {
        // The index of this object in the vector of AugmentedDeposit objects.
        int Index;
        // The number of photo electrons generated in the cube.
        double Deposit;
        // The attenuation between the cube and the MPPC.
        double Attenuation;
        // The cube in the vector of AugmentedCubes;
        std::weak_ptr<AugmentedCube> Cube;
        // The key of the fiber in the map of AugmentedFibers;
        std::weak_ptr<AugmentedFiber> Fiber;
        // The energy deposit going into this fiber from the associated cube.
        // This is not attenuation corrected.
        double GetDeposit() const {return Deposit;}
        // The energy deposit going into this fiber from the associated cube.
        // This is not attenuation corrected.
        void SetDeposit(double dep) {Deposit = dep;}
        // Change the energy deposit going into this fiber from the associated
        // cube. The change is in "pe".  This *does* *not* change the deposit
        // for other fibers contributing to the associated cube.  The change
        // is not attenuation corrected.
        void ChangeDeposit(double change);
        // Get contribution by this energy deposit to the measured charge of
        // this fiber.  This is attenuation corrected.
        double GetMeasurement() const {return Attenuation*Deposit;}
        // Set the contribution by this energy deposit to the expected
        // measured charge at the MPPC.
        void SetMeasurement(double measurement) {
            Deposit = measurement/Attenuation;
        }
        // A convenient way to check if the fiber is shared between deposits.
        bool HasSharedFiber() const;
        // A convenient way to get the charge actually measured by the fiber
        // associated with this deposit.
        double GetFiberMeasurement() const;
        // A convenient way to get the sum of expected measurements (from all
        // of the cubes) for the fiber.  This is the sum of all of the
        // expected measurements associated with the fiber.
        double GetFiberExpectedMeasurement() const;
        // The number of shared cubes on the fiber for this deposit.
        int GetFiberCubes() const;
        // Calculate the derivative for changing just *this* deposit.  This is
        // broken into two components.  The component for the fiber, and the
        // component for the cube.  Those two components are calculated in
        // separate methods.  The alpha parameter is a Lagrange multiplier.
        // The sum of the deposits on a fiber *must* add up to the
        // measurement, but that constraint makes the basic equations that are
        // being minimized singular.  The true minimum of the calculation will
        // meet the constraint.  The multiplier will start small, and be
        // increased until it's (approximately) infinite.
        double ConstraintDerivative(double alpha) const;
        // The contribution to the derivative from the fiber for changing
        // *just* this deposit in the the fiber.  The likelihood for the
        // measurement in the fiber is Poissonian, but this uses the Gaussian
        // approximation in all cases.
        double FiberDepositDerivative() const;
        // The contribution to the derivative from the cube for changing *just*
        // this deposit in the cube.
        double CubeDepositDerivative() const;
    };

    // Calculate the attenuation based on the distance the light is traveling
    // in a fiber.
    double Attenuation(Cube::Handle<Cube::Hit> hit, double dist);

    // Fill all of the augmented cubes, fibers and deposits.  The augmented
    // objects are used to track some extra information used by the charge
    // sharing calculation.  They are temporary, and all of the important
    // information is transfered to the CHit3D object at the end of the
    // calculation.
    void FillAugmented(const Cube::HitSelection& hit3D);

    // Do one "relaxation" step.  This is a way to minimize arbitrarily large
    // numbers of hits without breaking things.  It's slow, and basically an
    // implementation of steepest descent.  This is evolving the hits
    // according to the constraint that all of the fibers contributing to a
    // cube have the "same" contribution while the fiber charge is conserved
    // (with a lagrange multiplier).
    double EvolveConstraints(double step, double alpha);

    double EvolveCubes(double step, std::vector<double>& grad);

public:
    ShareCharge();
    ~ShareCharge();

    // This applies the charge sharing to the input THitSelection.  If you
    // don't want to change the original hits, then you must copy them first.
    void ApplyConstraints(Cube::HitSelection& mutableHits);

    void OptimizeCubes(Cube::HitSelection& mutableHits);

    // This applies charge sharing to the input THitSelection but minimizing
    // the chi2 against the fiber measurement.  For a given set of
    // constraints, it maximizes the charge entropy.  If you don't want to
    // change the original hits, then you must copy them first.
    void MaximizeEntropy(Cube::HitSelection& mutableHits);

    // Get the total deposit for all of the cubes.  This cheats by summing the
    // AugmentedDeposit values, and not accessing the cubes.
    double GetTotalDeposit();

    // Get the sum of the charge measured by the fibers.
    double GetTotalCharge();

    /// Get the sum of the expected charge for the fibers.
    double GetExpectedTotalCharge();

    // Get the total entropy for the cubes.  This is based on the idea that
    // all of the cubes should have the same deposit (the maximum entropy
    // criteria, similar to what is used in image processing).  It's adjusted
    // so that if every cube has the same deposit (i.e. the state with the
    // highest entropy), the total entropy is zero.  That means that the
    // totalEntropy is always less than or equal to zero.
    double GetTotalEntropy();

    // Get a fast estimate of the entropy for the cubes.  This cheats (big
    // time) by using a Gaussian approximation near the minimum entropy point.
    // This has the advantage that the gradiant is easily calcuable.
    double GetFastEntropy();

    // Calculate the partial derivative of the fast entropy approximation for
    // a single deposit.
    double GetFastEntropyDerivative(double d);

    // Get the likelihood for all of the fiber expected and measured charges.
    // This is the "ln(Poisson" likelihood, but using the Gaussian
    // approximation.
    double GetFiberChi2();

    // Get the number of augmented cubes.
    int GetAugmentedCubeCount() {return fAugmentedCubes.size();}

    // Get a reference to the augmented cube.  This is exposing the cubes so
    // that a non-member minimizer can work on stuff.
    AugmentedCube& GetAugmentedCube(int i);

    // Set a flag to indicate whether the reconstructed charge deposits in all
    // of the cubes should be rescaled to conserve the total measured charge
    // in the MPPCs.  If this is false, then the deposit will minimize the
    // chi2 for the charge distribution.  If this is true, then the deposit is
    // rescaled.  The rescaling will usually slightly increase the chi2 value
    // (typical rescaling is a few percent), but will make sure that the total
    // visible energy is not affected.
    void SetChargeConservation(bool v) {fChargeConservation = v;}

private:
    bool fChargeConservation;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
