#include "CubeShareCharge.hxx"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TRandom.h>

#include <set>

Cube::ShareCharge::ShareCharge() : fChargeConservation(true) {}
Cube::ShareCharge::~ShareCharge() {}

double Cube::ShareCharge::AugmentedCube::CubeDepositDerivative(
    double q1, double q2, double q3) const {
    // The derivative of the q1 deposit based on the other deposits in the
    // cube.  This is actually just a function, but is here since it is
    // conceptually associated with the cube.  If one of the fibers is
    // missing, then the associated charge should be set to zero.  The formula
    // was generated using wxmaxima, and is the derivative of
    //
    // Qavg : (q1+q2+q3)/3
    // X2: [(q1-Qavg)^2 + (q2-Qavg)^2 + (q3+Qavg)^2]/(Qavg + 1)
    // f90(diff(X2,q1));
    //
    // The code is then tweaked by hand to make sure that it fits C++.  The
    // function artificially imposes the constraint that all of the charges
    // are greater or equal to zero.
    //
    if (q1 < 0.0) q1 = 0.0;
    if (q2 < 0.0) q2 = 0.0;
    if (q3 < 0.0) q3 = 0.0;
    double deriv = (((-2.0)*(q3-(q3+q2+q1)/3.0))/3.0+((-2.0)*(q2-(q3+q2+q1)/3.0))/3.0+(4.0*(q1-(q3+q2+q1)/3.0))/3.0)/((q3+q2+q1)/3.0+1)-((std::pow(q3-(q3+q2+q1)/3.0,2.0)+std::pow(q2-(q3+q2+q1)/3.0,2.0)+std::pow(q1-(q3+q2+q1)/3.0,2.0))/std::pow((q3+q2+q1)/3.0+1.0,2.0))/3.0;
    // There are only two degrees of freedom (three for the charges minus one
    // for the average), but the derivative is calculated for three so that
    // all of the charges appear in the formula the same way.  The returned
    // derivative is reduced by 2/3 to fix the overcounting of the actual
    // number of degrees of freedom.
    return 2.0*deriv/3.0;
}

double Cube::ShareCharge::AugmentedCube::GetDeposit() const {
    double deposit = 0.0;
    for (std::vector<std::weak_ptr<AugmentedDeposit>>::const_iterator d
             = Deposits.begin();
         d != Deposits.end(); ++d) {
        deposit += d->lock()->GetDeposit();
    }
    return deposit;
}

void Cube::ShareCharge::AugmentedCube::SetDeposit(double dep) {
    double deposits = (double) Deposits.size();
    if (dep < 1.0) dep = 1.0;
    dep = dep / deposits;
    for (std::vector<std::weak_ptr<AugmentedDeposit>>::const_iterator d
             = Deposits.begin();
         d != Deposits.end(); ++d) {
        d->lock()->SetDeposit(dep);
    }
}

double Cube::ShareCharge::AugmentedFiber::GetExpectedMeasurement() const {
    double deposits = 0.0;
    for (std::vector<std::weak_ptr<AugmentedDeposit>>::const_iterator
             d = Deposits.begin();
         d != Deposits.end(); ++d) {
        deposits += d->lock()->GetMeasurement();
    }
    return deposits;
}

void Cube::ShareCharge::AugmentedDeposit::ChangeDeposit(double change) {
    while (std::abs(change) > 0.5*Deposit) change = 0.5*change;
    Deposit += change;
}

bool Cube::ShareCharge::AugmentedDeposit::HasSharedFiber() const {
    return Fiber.lock()->IsSharedFiber();
}

double Cube::ShareCharge::AugmentedDeposit::GetFiberMeasurement() const {
    return Fiber.lock()->GetMeasurement();
}

double Cube::ShareCharge::AugmentedDeposit::GetFiberExpectedMeasurement()
    const {
    return Fiber.lock()->GetExpectedMeasurement();
}

int Cube::ShareCharge::AugmentedDeposit::GetFiberCubes() const {
    return Fiber.lock()->Deposits.size();
}

double Cube::ShareCharge::AugmentedDeposit::ConstraintDerivative(double alpha)
    const {
    std::shared_ptr<AugmentedFiber> theFiber = Fiber.lock();
    if (!theFiber->IsSharedFiber()) {
        CUBE_ERROR << "Calculating derivative for a fiber with only one hit"
                   << std::endl;
        return 0.0;
    }
    double fiberDerivative = FiberDepositDerivative();
    double cubeDerivative = CubeDepositDerivative();
    double deriv = cubeDerivative
        + alpha*theFiber->GetMeasurement()*fiberDerivative;
    return deriv;
}

double Cube::ShareCharge::AugmentedDeposit::FiberDepositDerivative()
    const {
    std::shared_ptr<AugmentedFiber> theFiber = Fiber.lock();
    double expected = theFiber->GetExpectedMeasurement();
    double seen = theFiber->GetMeasurement();
    // This is actually the derivative of chi2.
    double deriv = 2.0*(expected-seen)/Attenuation/seen;
    // If one of the contributions is getting too small, then reduce it's
    // derivative.
    if (deriv > 0.0 && GetDeposit() < 3.0) {
        double delta = (3.0-GetDeposit());
        // delta(deposit==3) == 1 and delta(deposit<1) == 0;
        delta = 0.5*(2.0-delta);
        if (delta < 0.0) delta = 0.0;
        deriv *= std::sqrt(delta);
    }
    return deriv;
}

double Cube::ShareCharge::AugmentedDeposit::CubeDepositDerivative()
    const {
    std::shared_ptr<AugmentedCube> theCube = Cube.lock();
    std::vector<double> q;
    q.push_back(GetDeposit());
    for (std::vector<std::weak_ptr<AugmentedDeposit>>::const_iterator
             d = theCube->Deposits.begin();
         d != theCube->Deposits.end(); ++d) {
        if (d->lock()->Cube.lock() == theCube) continue;
        q.push_back(d->lock()->GetDeposit());
    }
    if (q.size() > 3) {
        CUBE_ERROR << "CUBES DO NOT HAVE FOUR FIBERS"
                   << std::endl;
        throw;
    }
    while (q.size() < 3) q.push_back(0.0);
    double deriv = theCube->CubeDepositDerivative(q[0],q[1],q[2]);
    return deriv;
}

double Cube::ShareCharge::Attenuation(
    Cube::Handle<Cube::Hit> fiber, double dist) {
    double f = fiber->GetProperty("Ratio12");
    double t1 = fiber->GetProperty("Atten1");
    double t2 = fiber->GetProperty("Atten2");
    double p0 = f * std::exp(-dist/t1) + (1.0-f)*std::exp(-dist/t2);
    return p0;
}


void Cube::ShareCharge::FillAugmented(const Cube::HitSelection& hits3D) {
    fAugmentedCubes.clear();
    fAugmentedFibers.clear();
    fAugmentedDeposits.clear();

    // Extract the simple hits from all of the input composite Hits.
    std::set<Cube::Handle<Cube::Hit>> fiberHits;
    for (Cube::HitSelection::const_iterator hit = hits3D.begin();
         hit != hits3D.end(); ++hit) {
        if ((*hit)->GetConstituentCount() < 3) {
            CUBE_ERROR << "3D cube hits need at least 3 contributors"
                       << std::endl;
            continue;
        }
        for (int i = 0; i<(*hit)->GetConstituentCount(); ++i) {
            Cube::Handle<Cube::Hit> fiber = (*hit)->GetConstituent(i);
            if (fiber->GetConstituentCount() > 0) {
                CUBE_ERROR << "Hit is not for a fiber" << std::endl;
                continue;
            }
            fiberHits.insert(fiber);
        }
    }

    // Create the augmented fibers.
    for (std::set<Cube::Handle<Cube::Hit>>::iterator fiber = fiberHits.begin();
         fiber != fiberHits.end(); ++fiber) {
        std::map<Cube::Handle<Cube::Hit>,
                 std::shared_ptr<AugmentedFiber>>::iterator
            fiberCheck = fAugmentedFibers.find((*fiber));
        if (fiberCheck != fAugmentedFibers.end()) {
            CUBE_ERROR << "Fiber already added to augmented fibers"
                       << std::endl;
            CUBE_LOG(0) << "   Id   "
                        << " " << fiberCheck->second->Hit->GetIdentifier()
                        << " " << (*fiber)->GetIdentifier() << std::endl;
            CUBE_LOG(0) << "   Time   "<<fiberCheck->second->Hit->GetTime()
                       << " " << (*fiber)->GetTime() << std::endl;
            CUBE_LOG(0) << "   Charge "<<fiberCheck->second->Hit->GetCharge()
                       << " " << (*fiber)->GetCharge() << std::endl;
            continue;
        }
        std::shared_ptr<AugmentedFiber> newFiber(new AugmentedFiber);
        newFiber->Hit = (*fiber);
        newFiber->Measurement = (*fiber)->GetCharge();
        fAugmentedFibers[newFiber->Hit] = newFiber;
    }

    // Fill the augmented cubes and augmented deposits.
    for (Cube::HitSelection::const_iterator cube = hits3D.begin();
         cube != hits3D.end(); ++cube) {
        std::shared_ptr<AugmentedCube> newCube(new AugmentedCube);
        newCube->Hit = *cube;
        for (int i = 0; i<(*cube)->GetConstituentCount(); ++i) {
            Cube::Handle<Cube::Hit> fiber = (*cube)->GetConstituent(i);
            double dist = ((*cube)->GetPosition() - fiber->GetPosition()).Mag();
            std::shared_ptr<AugmentedDeposit> newDeposit(new AugmentedDeposit);
            newDeposit->Cube = newCube;
            newDeposit->Cube.lock()->Deposits.push_back(newDeposit);
            newDeposit->Attenuation = Attenuation(fiber,dist);
            std::map<Cube::Handle<Cube::Hit>,
                     std::shared_ptr<AugmentedFiber>>::iterator
                fiberCheck = fAugmentedFibers.find(fiber);
            if (fiberCheck == fAugmentedFibers.end()) throw;
            newDeposit->Fiber = fiberCheck->second;
            newDeposit->Fiber.lock()->Deposits.push_back(newDeposit);
            newDeposit->SetMeasurement(newDeposit->Fiber.lock()->Measurement);
            fAugmentedDeposits.push_back(newDeposit);
        }
        fAugmentedCubes.push_back(newCube);
    }

    // Set the deposit measurements so that they sum to the total fiber charge.
    for (std::vector<std::shared_ptr<AugmentedDeposit>>::iterator d
             = fAugmentedDeposits.begin();
         d !=  fAugmentedDeposits.end(); ++d) {
        if ((*d)->GetFiberCubes()<1) continue;
        (*d)->SetMeasurement(
            (*d)->GetFiberMeasurement()/(*d)->GetFiberCubes());
    }

    // Check the number of fibers with overlaps.
    int overlaps = 0;
    for (std::map<Cube::Handle<Cube::Hit>,
             std::shared_ptr<AugmentedFiber>>::iterator f
             = fAugmentedFibers.begin();
         f !=  fAugmentedFibers.end(); ++f) {
        if (f->second->IsSharedFiber()) ++overlaps;
    }

    CUBE_LOG(0) << "Augmented Cubes      " << fAugmentedCubes.size() << std::endl;
    CUBE_LOG(0) << "Augmented Deposits   " << fAugmentedDeposits.size() << std::endl;
    CUBE_LOG(0) << "Augmented Fibers     " << fAugmentedFibers.size()
             << "   (" << overlaps << " with overlaps)" << std::endl;
    CUBE_LOG(0) << "Total Entropy        " << GetTotalEntropy() << std::endl;
    CUBE_LOG(0) << "Total Deposit        " << GetTotalDeposit() << std::endl;
    CUBE_LOG(0) << "Total Charge         " << GetTotalCharge() << std::endl;
}

double Cube::ShareCharge::GetTotalDeposit() {
    double totalDeposit = 0.0;
    for (std::vector<std::shared_ptr<AugmentedDeposit>>::iterator d
             = fAugmentedDeposits.begin();
         d !=  fAugmentedDeposits.end(); ++d) {
        totalDeposit += (*d)->GetDeposit();
    }
    return totalDeposit;
}

double Cube::ShareCharge::GetTotalCharge() {
    double totalCharge = 0.0;
    for (std::map<Cube::Handle<Cube::Hit>,
             std::shared_ptr<AugmentedFiber>>::iterator f
             = fAugmentedFibers.begin();
         f != fAugmentedFibers.end(); ++f) {
        totalCharge += f->second->GetMeasurement();
    }
    return totalCharge;
}

double Cube::ShareCharge::GetExpectedTotalCharge() {
    double totalCharge = 0.0;
    for (std::map<Cube::Handle<Cube::Hit>,
             std::shared_ptr<AugmentedFiber>>::iterator f
             = fAugmentedFibers.begin();
         f != fAugmentedFibers.end(); ++f) {
        totalCharge += f->second->GetExpectedMeasurement();
    }
    return totalCharge;
}

double Cube::ShareCharge::GetTotalEntropy() {
    double totalDeposit = GetTotalDeposit();
    double cubes = fAugmentedCubes.size();
    double totalEntropy = 0.0;
    if (cubes < 1) return 0.0;
    for (std::vector<std::shared_ptr<AugmentedCube>>::iterator c
             = fAugmentedCubes.begin();
         c !=  fAugmentedCubes.end(); ++c) {
        double deposit = (*c)->GetDeposit();
        double entropy = -deposit*std::log(deposit/totalDeposit)/totalDeposit;
        totalEntropy += entropy;
    }
    double avg = totalDeposit/cubes;
    // Shift so the maximum entropy is zero.
    double entropy = -cubes*avg*std::log(avg/totalDeposit)/totalDeposit;
    totalEntropy = totalEntropy - entropy;
    return totalEntropy;
    // Apropos nothing: For large number of cubes, the magnitude of the second
    // partial derivative with respect to a single cube deposit (when the
    // deposit is near the average deposit) is about
    // 1/totalDeposit/sqrt(cubes).
}

double Cube::ShareCharge::GetFastEntropy() {
    double totalDeposit = GetTotalDeposit();
    double cubes = fAugmentedCubes.size();
    double averageDeposit = totalDeposit/cubes;
    double fastEntropy = 0.0;
    if (cubes < 1) return 0.0;
    for (std::vector<std::shared_ptr<AugmentedCube>>::iterator c
             = fAugmentedCubes.begin();
         c !=  fAugmentedCubes.end(); ++c) {
        double deposit = (*c)->GetDeposit();
        double entropy = (deposit-averageDeposit);
        entropy = entropy*entropy;
        fastEntropy += entropy;
    }
    return - fastEntropy/std::sqrt(cubes)/totalDeposit;
}

double Cube::ShareCharge::GetFastEntropyDerivative(double d,
                                                   double totalDeposit) {
#ifdef SLOW_DOWN_PLEASE
    double total = GetTotalDeposit();
    if (abs(total-totalDeposit)>0.1) {
        CUBE_ERROR << "TOTAL  " << total << " " << totalDeposit << std::endl;
        throw std::runtime_error("Bad total deposit");
    }
#endif
    double cubes = fAugmentedCubes.size();
    double averageDeposit = totalDeposit/cubes;
    double deriv = - 2.0*(1.0-1.0/cubes)*(d-averageDeposit);
    deriv /= totalDeposit*std::sqrt(cubes);
    return deriv;
}

double Cube::ShareCharge::GetFiberChi2() {
    double chi2 = 0.0;
    for (std::map<Cube::Handle<Cube::Hit>,
             std::shared_ptr<AugmentedFiber>>::iterator f
             = fAugmentedFibers.begin();
         f != fAugmentedFibers.end(); ++f) {
        double m = f->second->GetMeasurement();
        double e = f->second->GetExpectedMeasurement();
        double d = (m-e);
        chi2 += d*d/m;
    }
    return chi2;
}

// Get a reference to the augmented cube.  This is exposing the cubes so
// that a non-member minimizer can work on stuff.
Cube::ShareCharge::AugmentedCube&
Cube::ShareCharge::GetAugmentedCube(int i) {
    if (i<0) throw std::runtime_error("Cube out of range");
    if (fAugmentedCubes.size() <= i) {
        throw std::runtime_error("Cube out of range");
    }
    return *fAugmentedCubes[i];
}

double Cube::ShareCharge::EvolveConstraints(double step, double alpha) {
    double totalChange = 0.0;
    int deposit = 0;
    for (std::vector<std::shared_ptr<AugmentedDeposit>>::iterator d
             = fAugmentedDeposits.begin();
         d !=  fAugmentedDeposits.end(); ++d) {
        // If this is the only deposit for the fiber, then the deposit is
        // fixed to the measurement on the fiber.
        if (!(*d)->HasSharedFiber()) {
            double r = (*d)->GetMeasurement() - (*d)->GetFiberMeasurement();
            if (std::abs(r) > 0.0001) {
                CUBE_ERROR << "Single fiber disagreement delta " << r
                           << " force " << (*d)->GetMeasurement()
                           << " to " << (*d)->GetFiberMeasurement()
                           << std::endl;
            }
            (*d)->SetMeasurement((*d)->GetFiberMeasurement());
            continue;
        }

        double deriv = (*d)->ConstraintDerivative(alpha);
        deriv = step * deriv / (*d)->GetFiberCubes();
        totalChange += std::abs(deriv);
#ifdef DEBUG_CHANGES
#undef DEBUG_CHANGES
        std::cout << "Constraint " << (*d)->GetDeposit()
                  << " contributes " << (*d)->GetMeasurement()
                  << " to " << (*d)->GetFiberMeasurement()
                  << " out of " << (*d)->GetFiberExpectedMeasurement()
                  << " from " << (*d)->GetFiberCubes()
                  << " change " << - deriv
                  << std::endl;
#endif
        // Never make a step of more than one photoelectron.
        const double maxChange = 1.0;
        if (deriv > maxChange) deriv = maxChange;
        if (deriv < -maxChange) deriv = -maxChange;
        const double minChange = 0.01;
        if (std::abs(deriv) < minChange) {
            // Never take a step of less than 0.01 pe.
            if (0.0 < deriv) deriv = minChange;
            else deriv = -minChange;
        }
        // The derivative points away from the minimum so "step backwards"
        (*d)->ChangeDeposit(-deriv);
    }
    return totalChange;
}

void Cube::ShareCharge::ApplyConstraints(Cube::HitSelection& mutableHits) {
    FillAugmented(mutableHits);

    // Relax for a "very long time".  This could be a lot more efficient, but
    // it's not so slow, so WTH.
    double alpha = 0.1;
    int smallChi2 = 0;
    for (int i=0; i<10000; ++i) {
        double step = 0.5/alpha;
        // Take one step.
        double totalChange = EvolveConstraints(step,alpha);
        // Find the total measured charge, and the current sum of the
        // charge distributed to the cubes.  When the measured charge and
        // the distributed charge are close, stop the iterations.
        double measuredCharge = 0.0;
        double fiberCharge = 0.0;
        for (std::vector<std::shared_ptr<AugmentedDeposit>>::const_iterator
                 d = fAugmentedDeposits.begin();
             d != fAugmentedDeposits.end(); ++d) {
            measuredCharge += (*d)->GetMeasurement();
            fiberCharge += (*d)->GetFiberMeasurement()/(*d)->GetFiberCubes();
        }
        double diff = std::abs(measuredCharge - fiberCharge);
        double delta = diff/fiberCharge;
        double chi2 = GetFiberChi2();
#ifdef DEBUG_EVOLUTION
#undef DEBUG_EVOLUTION
        CUBE_LOG(0) << i
                 << " " << totalChange
                 << " " << alpha
                 << " " << step
                 << " " << delta
                 << " " << GetTotalEntropy()
                 << " " << GetFiberChi2() << std::endl;
#endif
        ++smallChi2;
        if (chi2 > 0.01) smallChi2 = 0.0;
        if (smallChi2 > 20) break;
        if (diff < 1.0 && delta < 1.0E-4) break;
        // Increase the Lagrange multiplier.
        alpha = std::min(1.001*alpha,100.0);
    };

    CUBE_LOG(0) << "Constrainted Deposit: " << GetTotalDeposit()
             << " Expected Q: " << GetExpectedTotalCharge()
             << " Orig: " << GetTotalCharge()
             << " Chi2 =" << GetFiberChi2() << std::endl;

    if (fChargeConservation) {
        double chargeRatio = GetTotalCharge()/GetExpectedTotalCharge();
        for (int i=0; i<GetAugmentedCubeCount(); ++i) {
            double q = GetAugmentedCube(i).GetDeposit();
            GetAugmentedCube(i).SetDeposit(chargeRatio*q);
        }
        CUBE_LOG(0) << "   Rescale Deposit: " << GetTotalDeposit()
                 << " Scaled by " << chargeRatio
                 << " Chi2 =" << GetFiberChi2()
                 << " S=" << GetTotalEntropy() << std::endl;
    }

    // Save the values into the reconstructed hits.
    for (std::vector<std::shared_ptr<AugmentedCube>>::iterator
             c = fAugmentedCubes.begin();
         c != fAugmentedCubes.end(); ++c) {
        Cube::Handle<Cube::WritableHit> hit = (*c)->Hit;
        if (!hit) {
            CUBE_ERROR << "The input hits need to be Writable<blah> objects"
                       << std::endl;
            continue;
        }
        hit->SetCharge((*c)->GetDeposit());
    }

}


double Cube::ShareCharge::EvolveCubes(double step,
                                      std::vector<double>& grad) {
    double totalDeposit = GetTotalDeposit();
    double averageDeposit = totalDeposit/fAugmentedCubes.size();
    if (grad.size () != fAugmentedCubes.size()) {
        grad.resize(fAugmentedCubes.size());
    }
    // Find the gradient.
    double magGrad = 0.0;
    for (int c = 0; c < GetAugmentedCubeCount(); ++c) {
        AugmentedCube& cube = GetAugmentedCube(c);
        double currentDeposit = cube.GetDeposit();
        double fiberDeposit = currentDeposit/cube.Deposits.size();
        cube.SetDeposit(currentDeposit); // force all deposits to be the same.
        // Find the derivative of the deposit. This is the sum of the
        // derivatives for all fibers.
        double fiberDeriv = 0;
        for (std::vector<std::weak_ptr<AugmentedDeposit>>::iterator d
                 = cube.Deposits.begin();
             d != cube.Deposits.end(); ++d) {
            double dir = d->lock()->FiberDepositDerivative();
            fiberDeriv += dir;
        }
        // Find the derivative due to the change in the entropy (and apply a
        // scale so that it doesn't get lost in the noise.  This might be
        // needed due to a math error, or might just be needed.
        double eScale = 1.0*fAugmentedCubes.size();
        double entropyDeriv = eScale*GetFastEntropyDerivative(currentDeposit,
                                                              totalDeposit);
        double deriv = fiberDeriv - entropyDeriv;
        grad[c] = deriv;
        magGrad += deriv*deriv;
    }
    magGrad = std::max(std::sqrt(magGrad),0.00001);
    // Make the step.
    double dist = 0.0;
    for (int c = 0; c < GetAugmentedCubeCount(); ++c) {
        AugmentedCube& cube = GetAugmentedCube(c);
        // The derivative points away from the minimum so "step backwards"
        double change = - step * grad[c]/magGrad;
        double currentDeposit = cube.GetDeposit();
        double targetDeposit = currentDeposit + change;
        cube.SetDeposit(targetDeposit);
#ifdef DEBUG_CHANGES
#undef DEBUG_CHANGES
        double newDeposit = cube.GetDeposit();
        if (std::abs(targetDeposit - newDeposit) > 1E-4) {
            continue;
        }
        dist += change*change;
        std::cout << "Optim Cube " << cube.GetDeposit()
                  << " " << grad[c]
                  << " " << change
                  << std::endl;
#endif
    }
    double chi2 = GetFiberChi2();
    double entropy = GetFastEntropy();
    return chi2 + entropy;
}

void Cube::ShareCharge::OptimizeCubes(Cube::HitSelection& mutableHits) {
    FillAugmented(mutableHits);

    // Relax for a "very long time".  This could be a lot more efficient, but
    // it's not so slow, so WTH.
    double step = std::sqrt(GetTotalDeposit()/GetAugmentedCubeCount());
    std::vector<double> grad;
    std::vector<double> oldGrad;
    double oldChi2 = -1.0;
    double oldChange = 0.0;
    int stuck = 0;
    for (int i=0; i<1000; ++i) {
        // Take one step.
        double chi2 = EvolveCubes(step,grad);
        double change = 0.0;
        double dotProd = 0.0;
        if (grad.size() != oldGrad.size()) {
            oldGrad.resize(grad.size());
            for (int i=0; i<grad.size(); ++i) {
                dotProd += grad[i]*grad[i];
                change += grad[i]*grad[i];
            }
            oldChange = change;
        }
        else {
            for (int i=0; i<grad.size(); ++i) {
                dotProd += grad[i]*oldGrad[i];
                change += grad[i]*grad[i];
            }
        }
        double cosGrad = dotProd/sqrt(change)/sqrt(oldChange);
        double deltaChi2 = chi2-oldChi2;
        ++stuck;
        if (oldChi2 < 0) deltaChi2 = -deltaChi2;
        else if (std::abs(deltaChi2) > 0.01) stuck = 0;
        oldChi2 = chi2;
        std::copy(grad.begin(),grad.end(),oldGrad.begin());
        oldChange = change;
#ifdef DEBUG_EVOLUTION
#undef DEBUG_EVOLUTION
        CUBE_LOG(0) << i
                    << " f " << chi2
                    << " (" << deltaChi2 << "," << stuck << ")"
                    << " s " << step
                    << " g " << dotProd << " " << cosGrad << std::endl;
#endif
        // This approximates a golden-section search for the minimum...
        if (stuck > 5)  break;
        else if (deltaChi2 > 0) step = step/1.618;
        else if (step < 0.05) break;
        else if (dotProd > 0.0) step = 1.15*step;
    };

    CUBE_LOG(0) << "Optimized Deposit: " << GetTotalDeposit()
             << " Expected Q: " << GetExpectedTotalCharge()
             << " Orig: " << GetTotalCharge()
             << " Chi2 =" << GetFiberChi2() << std::endl;

    if (fChargeConservation) {
        double chargeRatio = GetTotalCharge()/GetExpectedTotalCharge();
        for (int i=0; i<GetAugmentedCubeCount(); ++i) {
            double q = GetAugmentedCube(i).GetDeposit();
            GetAugmentedCube(i).SetDeposit(chargeRatio*q);
        }
        CUBE_LOG(0) << "   Rescale Deposit: " << GetTotalDeposit()
                 << " Scaled by " << chargeRatio
                 << " Chi2 =" << GetFiberChi2()
                 << " S=" << GetTotalEntropy() << std::endl;
    }

    // Save the values into the reconstructed hits.
    for (std::vector<std::shared_ptr<AugmentedCube>>::iterator
             c = fAugmentedCubes.begin();
         c != fAugmentedCubes.end(); ++c) {
        Cube::Handle<Cube::WritableHit> hit = (*c)->Hit;
        if (!hit) {
            CUBE_ERROR << "The input hits need to be Writable<blah> objects"
                       << std::endl;
            continue;
        }
        hit->SetCharge((*c)->GetDeposit());
    }
}

namespace {
    // This fits the energy deposit in each cube to minimize the chi2 of the
    // fiber charge measurements, with the added criteria that the entropy for
    // the distribution of cube charges should be maximized.  The entropy
    // criteria breaks the degeneracy when there isn't a solution.  It imposes
    // a condition that (as much as possible) the cube deposits should all be
    // equal.  The input parameters are given in the log of the deposit so
    // that they are constrainted to be positive.  It's constructed so that
    // the ideal starting point is basically par[0..i] = 0.0;
    class CubeMaxEntropy : public ROOT::Math::IMultiGenFunction {
        Cube::ShareCharge* fShareCharge;
        double fAverageDeposit;
    public:
        explicit CubeMaxEntropy(Cube::ShareCharge* shared)
            : fShareCharge(shared) {
            fAverageDeposit = fShareCharge->GetTotalDeposit();
            fAverageDeposit /= fShareCharge->GetAugmentedCubeCount();
            CUBE_LOG(0) << "Average deposit is " << fAverageDeposit << std::endl;
        }
        ~CubeMaxEntropy() {}
        IBaseFunctionMultiDimTempl* Clone() const {return NULL;}
        unsigned int NDim() const {
            return fShareCharge->GetAugmentedCubeCount();
        }
        double DoEval(const double* par) const {
            for (int i=0; i<fShareCharge->GetAugmentedCubeCount(); ++i) {
                double v = fAverageDeposit * std::exp(par[i]);
                fShareCharge->GetAugmentedCube(i).SetDeposit(v);
            }
            double chi2 = fShareCharge->GetFiberChi2();
            double entropy = fShareCharge->GetTotalEntropy();
            double entropyScale = fShareCharge->GetAugmentedCubeCount();
            double goodness = chi2 - entropyScale*entropy;

            return goodness;
        }
    };
}

void Cube::ShareCharge::MaximizeEntropy(Cube::HitSelection& mutableHits) {
    FillAugmented(mutableHits);

    // Protect against really small events.
    if (GetAugmentedCubeCount() < 4) return;

    // Save the original values.
    double pre[10000];
    for (int i=0; i< GetAugmentedCubeCount(); ++i) {
        pre[i] = GetAugmentedCube(i).GetDeposit();
    }

    // Minimize!
    std::unique_ptr<ROOT::Math::Minimizer>
        minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    CubeMaxEntropy entropy(this);
    minimizer->SetFunction(entropy);
    for (int i=0; i<GetAugmentedCubeCount(); ++i) {
        std::ostringstream nm;
        nm << "Cube[" << i << "]";
        minimizer->SetVariable(i,nm.str(),0.0,0.10);
    }

    minimizer->SetPrintLevel(1);
    minimizer->Minimize();

    // Reevaluate at the standard deviation plus one.
    double unc[10000];
    for (int i=0; i< GetAugmentedCubeCount(); ++i) {
        unc[i] = minimizer->X()[i] + minimizer->Errors()[i];
    }
    entropy.DoEval(unc);
    for (int i=0; i< GetAugmentedCubeCount(); ++i) {
        unc[i] = GetAugmentedCube(i).GetDeposit();
    }

    double minChi2 = entropy.DoEval(minimizer->X());

    CUBE_LOG(0) << "MaxEnt Deposit: " << GetTotalDeposit()
             << " Total Q: " << GetTotalCharge()
             << " Expected Q: " << GetExpectedTotalCharge()
             << " Chi2 =" << GetFiberChi2()
             << " S=" << GetTotalEntropy()
             << " Minimum Value=" << minChi2 << std::endl;

    if (fChargeConservation) {
        double chargeRatio = GetTotalCharge()/GetExpectedTotalCharge();
        for (int i=0; i<GetAugmentedCubeCount(); ++i) {
            double q = GetAugmentedCube(i).GetDeposit();
            GetAugmentedCube(i).SetDeposit(chargeRatio*q);
        }
        CUBE_LOG(0) << "   Rescale Deposit: " << GetTotalDeposit()
                 << " Scaled by " << chargeRatio
                 << " Chi2 =" << GetFiberChi2()
                 << " S=" << GetTotalEntropy() << std::endl;
    }

    // Make sure the last evaluation is at the minimum, and save the values
    // into the reconstructed hits.
    for (int i=0; i< GetAugmentedCubeCount(); ++i) {
        Cube::Handle<Cube::WritableHit> hit = GetAugmentedCube(i).Hit;
        if (!hit) {
            CUBE_ERROR << "The input hits need to be Writable<blah> objects"
                       << std::endl;
            continue;
        }
        hit->SetCharge(GetAugmentedCube(i).GetDeposit());
        // hit->SetChargeUncertainty(unc[i]);
    }

}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
