#ifndef __SIMPLE_SIR_HH_SEEN__
#define __SIMPLE_SIR_HH_SEEN__

#include <TRandom.h>

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

// Define SIMPLE_SIR_NAMESPACE before including this code to place all of this
// into a namespace.
#ifdef SIMPLE_SIR_NAMESPACE
namespace SIMPLE_SIR_NAMESPACE {
#endif

/////////////////////////////////////////////////////////////////
// A simple template implementation of a Sequential Importance Resampling
// (SIR) Particle Filter.  I wrote this as an exercise to understand SIR
// particle filters, but it turns out to be surprisingly useful.  The file is
// almost entirely documentation since the actual template is only about 70
// lines of code.
//
// This defines a template which generates a class to apply a SIR particle
// filter.  The template operates on a weighted point cloud that describes a
// PDF.  The point cloud is defined as a vector of weight and state pairs (the
// StatePair class).  The user must construct a point cloud describing the
// prior information before the filter is applied.  Then the UpdateSamples()
// method is called for each new measurement.  When UpdateSamples() returns,
// the point cloud has been updated to reflect the new information.  See the
// example below.  There is also a working example included as an "ifdef"
// section of code at the end of this file.
//
// The template defines one important class method:
//
// SimpleSIR<>::UpdateSamples(SampleVector& samples, Measurement measurement)
// -- Apply the effects of the measurement to the PDF described by samples.
// The vector of samples is modified.
//
// The template defines classes that are useful for the user.
//
// SimpleSir<>::SamplePair -- A std::pair<double,UserState> that holds the
// weight (in first), and the user defined state in second.
//
// SimpleSir<>::SampleVector -- A std::vector<SamplePair> that holds the point
// cloud.
//
// The template takes four arguments.
//
// UserState -- A class that describes the state.  The particle filter doesn't
// try to access the internals of the state, so it can contain almost
// anything.  However, it must be copyable since the filter will copy the
// state by value.
//
// UserStatePropagator -- A class that implements ```double operator ()``` to
// propagate the state value to the time step (usually determined by the next
// measurement).  The return value should (almost) always be 1.0.  See below
// for details on how it needs to be declared.
//
// UserMeasurement -- A class that describes a measurement.  The particle
// filter doesn't try to access the internals of the measurement, so it can
// contain almost anything.  The template always uses the measurement as
// "const UserMeasurement&"
//
// UserLikelihood -- A class that implements ``` double operator ()`` to
// calculate the likelihood of a measurement given a state.  The return value
// needs to be between 0.0 and 1.0.  See below for details on how it needs to
// be declared.
//
// ## THE UserStatePropagator CLASS
//
// This is provided by the user to update a particular state to the next step.
// As far as the template is concerned, the class only needs to provide a
// single method declared.
//
// ```
// class UserStatePropagator {
// public:
//    double operator () (UserState& state,
//                        const UserMeasurement& measurement);
// };
// ```
//
// The return value is a weight adjustment that can be applied to the sample
// weight.  The return value should almost always be 1.0.  It is there so that
// the propagated state can have extra variance so that the propagated state
// does not get stuck in a "low variance" situation.  If you find that it's a
// problem, you probably should not use a particle filter (so always return
// 1.0).
//
// See below for an example implementation.
//
// ## The UserLikelihood CLASS
//
// This is provided by the user to calculate the likelihood of a measurement
// given a state.  As far as the template is concerned, the class only needs
// to provide a single method
//
// ```
// class UserLikelihood {
// public:
//     double operator() (const UserState& state,
//                        const UserMeasurement& measurement);
// };
// ```
//
// The return value needs to be between 0.0 and 1.0.
//
// ## EXAMPLE
//
// The following example shows how to implement a filter for measurements X_i
// made at T_i with a Gaussian uncertainty on the measurement of 1.0.  The
// true relation is "X_i ~ 1.0*T_i + E" where E is an error term due to noise.
//
// The X position is measured at a particular time t.  It's declared as a
// struct so that we don't need to add "public".  The default copy will be
// good enough.
//
// ```
// struct MyMeasure {
//     float t;
//     float x;
// };
// ```
//
// The state is the current position at a particular time.
//
// ```
// struct MyState {
//     float x;
//     float t;
// };
// ```
//
// The measurement has a unit uncertainty, so this just returns the
// unnormalized likelihood.  The weights are going to be normalized, so the
// likelihood doesn't need to be.
//
// ```
// struct MyLikelihood {
//    double operator () (const MyState& s, const MyMeasure& m) {
//        double d = m.x - s.x;
//        return std::exp(d*/2.0);
//    }
// };
// ```
//
// The propagator needs to advance the state to the new time, and then add the
// effect of noise to the propagation.  This should always return 1.0!
//
// ```
// struct MyPropagator {
//     double operator () (MyState& s,
//                         const MyMeasure& m) {
//         s.x = s.x + 1.0*(m.t - s.t);        /* Update the state */
//         s.x = s.x + gRandom->Gaus(0.0,0.1); /* Add the noise */
//         s.t = m.t;
//         return 1.0;
//     }
// };
// ```
//
// Declare a typedef for the template!
//
// ```
// typedef SimpleSIR<MyState, MyPropagator,
//                   MyMeasurement, MyLikelihood> MySIR;
// ```
//
// Assume that you've got a vector of measurements from someplace
//
// ```
// std::vector<MyMeasure> MeasurementVector;
// ```
//
// You will need to make a vector of states that describes your prior
// knowledge.
//
// ```
// double t0 = MeasurementVector.front().t;
// double x0 = MeasurementVector.front().x;
// MySIR::SampleVector samples;
// for (int i=0; i<1000; ++i) {
//    MySIR::SamplePair sample;
//    sample.first = 1.0;    /* "Always" 1.0 */
//    sample.second.x = gRandom(x0,1.0);
//    sample.second.t = t0;
//    samples.push_back(sample);
// };
// ```
//
// Now update through all of the measurements.  Notice that the first
// measurement is skipped since I used it in the prior.
//
// ```
// MySIR mySIR;
// for (std::size_t i = 1; i<MeasurementVector.size(); ++i) {
//    mySIR.UpdateSamples(samples,MeasurementVector[i]);
// }
// ```
//
// After applying all of the measurements, the final estimator for the state
// can be calculated as
//
// ```
// double avgX = 0.0;
// for (std::size_t i = 0; i<samples.size(); ++i) {
//     avgX += samples[i].first * samples[i].second.x;
// }
// ```
template<typename UserState, typename UserStatePropagator,
         typename UserMeasurement, typename UserLikelihood>
class SimpleSIR {
public:
    // Add a typedef for the UserState.  It makes the template code a little
    // easier, but external user code probably wants to use the real class
    // name.
    typedef UserState State;

    // Add a typedef for the UserMeasurement.  It makes the template code a
    // little easier, but external user code probably wants to use the real
    // class name.
    typedef UserMeasurement Measurement;

    // Create a field holding the class to propagate the state forward in
    // "time".  The class must provide a public method:
    //
    // ```C++
    // double operator () (UserState& s, const UserMeasurement& m);
    // ```
    //
    // The state, s, is modified to correspond to the measurement, m,
    // (i.e. moved to the closest position or updated to the time of the
    // measurement).  The return value is a weight adjustment to be applied to
    // the state (if you don't understand what this is for, return 1.0).  See
    // the example in the main documentation.
    UserStatePropagator Propagator;

    // Create a field holding the class to calculate the likelihood.  The
    // class must provide a public method:
    //
    // ```C++
    // double operator() (const UserState& state,
    //                    const UserMeasurement& measurement);
    // ```
    //
    // The return value is the likelihood of the measurement given the state,
    // and must be between 0.0 and 1.0.
    UserLikelihood Likelihood;

    // Each sample is a pair that consists of the weight for the state, and
    // the state.  The first element of the pair is the weight.  The second
    // element is the state.
    typedef std::pair<double,State> SamplePair;

    // A vector of all of the samples that represent the PDF for the state.
    // This is a cloud of states where the density and weights represent the
    // PDF, and the user will need to provide a vector of this type to the
    // UpdateSamples method.  The mean of the PDF is the sum of
    // SamplePair::first times SamplePair::second.  See the main comment above
    // for an example.
    typedef std::vector<SamplePair> SampleVector;

    // An iterator through the sample.  User code probably wants to directly
    // use "SampleVector::iterator", but this prevents the template from
    // needed typename every time the iterator is used.
    typedef typename SampleVector::iterator SampleIterator;

    // The constructor will set some default values.  It doesn't do much.
    SimpleSIR() {SetResampleFraction(0.5);}

    // Update the sample vector of states based on the provided measurement.
    // This will be called once by the user as each new measurement is added
    // to the fit.  Since the vector of states is modified, the user should
    // record any information before moving to the next measurement.  That
    // usually means that the user will want to record the sample mean and
    // variance, or possibly even copy the states to another vector.
    //
    // This returns false if there was no resampling, and true if there was a
    // resampling.  This will be called by the user for each step
    bool UpdateSamples(SampleVector& samples,
                       const Measurement& measurement) {

        // PROPAGATE: This propagates each sample to the measurement.  It
        // usually adds "noise" or multiple scattering.  The propagator should
        // return a weight adjustment (almost always 1.0).  This can be used
        // if the propagator is applying some excess variability so that the
        // states don't get "stuck".  The weight adjustment may be larger than
        // 1.0.
        for (SampleIterator s = samples.begin();
             s != samples.end(); ++s) {
            s->first *= Propagator(s->second, measurement);
        }

        // MEASURE: The finds the likelihood for each state given the
        // measurement.
        double norm = 0.0;
        for (SampleIterator s = samples.begin();
             s != samples.end(); ++s) {
            s->first *= Likelihood(s->second, measurement);
            norm += s->first;  // Find the normalization sum for later.
        }

        // RENORMALIZE: Everything to sums to 1
        for (SampleIterator s = samples.begin(); s != samples.end(); ++s) {
            s->first /= norm;
        }

        // EFFECTIVE SAMPLE SIZE: This checks if effective number of samples
        // is still large enough, and then might resample.  The effective
        // number of samples is 1.0/(variance of the weights).
        double effective = GetEffectiveSamples(samples.begin(), samples.end());

        // RESAMPLE: Resample if the effective sample size is too small.
        if (effective > ResampleFraction*samples.size()) return false;

        // If necessary, resize the work area.  This should only happen once.
        if (ResampleWorkArea.size() != samples.size()) {
            ResampleWorkArea.resize(samples.size());
        }

        // Sort the samples by weight.  This makes the resampler more
        // efficient since it checks the high weight samples first.
        std::sort(samples.begin(),samples.end());

        // And apply the resampling.
        for (SampleIterator s = ResampleWorkArea.begin();
             s != ResampleWorkArea.end(); ++s) {
            *s = Resample(samples.begin(),samples.end());
            s->first = 1.0/samples.size(); // uniform weight after resampling
        }

        // Replace with the new sample.
        samples = ResampleWorkArea;

        return true;
    }

    // The effective number of samples gets reduced as the filter runs.  When
    // the ratio of effective samples to actual samples drops below this
    // fraction, the filter will resample.  Set this to 1.0 to resample every
    // time, or 0.0 to never resample.
    void SetResampleFraction(double r) {
        if (r > 1.0) ResampleFraction = 1.0;
        if (r < 0.0) ResampleFraction = 0.0;
        ResampleFraction = r;
    }

    // Calculate the effective number of samples. The effective number of
    // samples is 1.0/(variance of the weights).
    double GetEffectiveSamples(SampleIterator begin, SampleIterator end) {
        double weightVariance = 0.0;
        while (begin != end) {
            weightVariance += begin->first * begin->first;
            ++begin;
        }
        return 1.0 / weightVariance;
    }

private:
    // When the effective number of particles has dropped below this fraction,
    // of the real number of samples, then resample.
    double ResampleFraction;

    // Pull one sample from the sample vector based on the sample weights.
    // The input samples will be sorted by weight so that the high weight
    // samples are at the end.  If this is to slow, it could be implemented
    // with a binary search (but the weights would need to be summed first).
    SamplePair Resample(SampleIterator begin,
                        SampleIterator end) {
        double norm = gRandom->Uniform(0.0,1.0);
        while (begin < end && norm > 0) {
            --end;
            norm -= end->first;
        }
        return *end;
    }

    // Provide a work area vector for resampling.  This prevents it from being
    // constructed and deconstructed on each resampling.  It will be resized
    // (hopefully, only once) during the first resampling.
    SampleVector ResampleWorkArea;

};

// MIT License

// Copyright (c) 2017-2020 Clark McGrew
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#ifdef SIMPLE_SIR_NAMESPACE
};
#endif

///////////////////////////////////////////////////////////////////////////
// End of SimpleSIR
///////////////////////////////////////////////////////////////////////////
#endif

#ifdef SIMPLE_SIR_HH_ADD_WORKING_EXAMPLE
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// WORKING EXAMPLE OF "TRACKFIT" WITH POSITION MEASUREMENTS
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// This next section is a more complicated working example used to debug this
// code using ROOT6.  It can be be run as test code in root 6 by doing
//
// cat RunSimpleSIRExample.C << EOF
// #define SIMPLE_SIR_HH_ADD_WORKING_EXAMPLE
// #include "SimpleSIR.hh"
// EOF
//
// and then running
//
// root RunSimpleSIRExample.C++
//
// It produces RunSimpleSIRExample.png.
//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include <TRandom.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>

// State 0: x, 1: y, 2: z, 3: dx, 4: dy, 5: dz
typedef std::vector<float> ExampleState;

// Measurement 0: x, 1: y, 2: z
typedef std::vector<float> ExampleMeasurement;

// Calculate the likelihood for the measurement given the state.  Provided by
// the user.
struct ExampleLikelihood {
    double operator() (const ExampleState& state,
                       const ExampleMeasurement& measurement) {
#ifdef NON_GAUSSIAN_LIKELIHOOD
        double maxpos = 0.0;
        for (int i=0; i<3; ++i) {
            double r = state[i] - measurement[i];
            if (r < 0.0) r = -r;
            if (maxpos < r) maxpos = r;
        }
        if (maxpos < 0.5) return 0.99;
        return 0.01;
#endif
        // GAUSSIAN (A fallback for the test).
        double r1 = 0.0;
        for(int i=0; i<3; ++i) {
            double r = state[i] - measurement[i];
            r1 += r*r;
        }
        return std::exp(-r1/2.0);
    }
};

// Update the state to a new measurement.  The measurement can be ignored, but
// is provided in case it's needed.  Provided by the user.  This example
// assumes 0) x, 1) y, 2) z, 3) dx, 4) dy, 5) dz
struct ExamplePropagate {
    double operator () (ExampleState& state,
                        const ExampleMeasurement& measurement) {
        // Make an initial estimate of the distance to closest approach.
        double dist = 0.0;
        for (int i=0; i<3; ++i) {
            dist += state[i+3]*(measurement[i]-state[i]);
        }
        // Add scattering to the direction and then normalize.
        double norm = 0.0;
        for (int i=3; i<6; ++i) {
            state[i] += dist*gRandom->Gaus(0.0,0.05);
            norm += state[i]*state[i];
        }
        norm = std::sqrt(norm);
        for (int i=3; i<6; ++i) {
            state[i] /= norm;
        }
        // Update the estimate of the distance to closest approach.
        dist = 0.0;
        for (int i=0; i<3; ++i) {
            dist += state[i+3]*(measurement[i]-state[i]);
        }
        // Update the position.
        for (int i=0; i<3; ++i) {
            state[i] += dist*state[i+3] + dist*gRandom->Gaus(0.0,0.1);
        }
        return 1.0;
    }
};

// Declare the filter.
typedef SimpleSIR<ExampleState, ExamplePropagate,
                  ExampleMeasurement, ExampleLikelihood> ExampleSIR;

// Calculate the average and covariance for the current state vector.  This is
// convenience function the user might want to write.
void ExampleMakeAverage(const ExampleSIR::SampleVector& samples,
                        ExampleState& stateAverage,
                        TMatrixD& stateCov) {
    std::size_t dim = samples[0].second.size();
    if (stateAverage.size() != dim) {
        std::cout << "Resize!!!" << std::endl;
        stateAverage.resize(dim);
        stateCov.ResizeTo(dim,dim);
    }
    for (std::size_t i=0; i<dim; ++i) {
        stateAverage[i] = 0.0;
        for (std::size_t j=0; j<dim; ++j) {
            stateCov(i,j) = 0.0;
        }
    }
    // Find the averages.
    for (ExampleSIR::SampleVector::const_iterator s = samples.begin();
             s != samples.end(); ++s) {
        for (std::size_t i=0; i<dim; ++i) {
            stateAverage[i] += s->first*s->second[i];
        }
    }
    // Find the covariance (brute force!)
    for (ExampleSIR::SampleVector::const_iterator s = samples.begin();
             s != samples.end(); ++s) {
        for (std::size_t i=0; i<dim; ++i) {
            double a = s->second[i]-stateAverage[i];
            for (std::size_t j=0; j<dim; ++j) {
                double b = s->second[j]-stateAverage[j];
                stateCov(i,j) += s->first*a*b;
            }
        }
    }
}

void RunSimpleSIRExample() {
    int sampSize = 10000;
    TH2F* hist1 = new TH2F("hist1","The PDF for the first two coordinates",
                           100, -5.0, 5.0,
                           100, -5.0, 5.0);

    // Fill a vector of measurements.  This needs to be provided by the user.
    std::vector<ExampleMeasurement> measurements;
    for (int i = 0; i < 100; ++i) {
        ExampleMeasurement m(4);
#ifdef RANDOM_MEASUREMENTS
        m[0] = gRandom->Gaus(); // x
        m[1] = gRandom->Gaus(); // y
#else
        m[0] = 0.0; // x
        m[1] = 0.0; // y
#endif
        m[2] = 1.0*i; // z
        m[3] = 1.0*i; // t
        measurements.push_back(m);
    }

    // Fill the sample vector with a prior distribution.  This needs to be
    // provided by the user.
    ExampleSIR::SampleVector samples;
    samples.clear();
    for (int i=0; i < sampSize; ++i) {
        ExampleState s(6);
        s[0] = measurements[0][0] + gRandom->Gaus(); // prior for x
        s[1] = measurements[0][1] + gRandom->Gaus(); // prior for y
        s[2] = measurements[0][2];                   // prior for z
        s[3] = measurements[1][0] + gRandom->Gaus(); // next x
        s[4] = measurements[1][1] + gRandom->Gaus(); // next y
        s[5] = measurements[1][2];                   // next z
        for (int j = 0; j<3; ++j) s[j+3] = s[j+3] - s[j];
        double r = 0.0;
        for (int j = 3; j<6; ++j) r += s[j]*s[j];
        for (int j = 3; j<6; ++j) s[j] /= r;
        ExampleSIR::SamplePair sample(1.0,s);
        samples.push_back(sample);
    }

#ifdef AllowInitialImportanceSampling
    // Correct for biased importance sampling.  Only needed if the initial
    // state weights are not all 1.0.  This almost certainly shouldn't be
    // done, but if you need it, here it is.  Don't do this unless you KNOW
    // what your are doing.
    for (ExampleSIR::SampleIterator s = samples.begin();
         s != samples.end(); ++s) {
        s->first = 1.0/s->first;
    }
#endif


    ///////////////////////////////////////////
    // Filter code starts here.
    ///////////////////////////////////////////

    // Declare the filter.
    ExampleSIR Filter;
    // Update the PDF described by the samples to the next measurement.
    ExampleState stateAverage;
    TMatrixD stateCovariance;
    for (std::size_t i = 1; i < measurements.size(); ++i) {
        Filter.UpdateSamples(samples,measurements[i]);
        // Calculate and print the averages.
        ExampleMakeAverage(samples,stateAverage,stateCovariance);
        for (ExampleState::iterator s = stateAverage.begin();
             s != stateAverage.end(); ++s) {
            std::cout << *s << " ";
        }
        std::cout << std::endl;
        stateCovariance.Print();
    }

    ///////////////////////////////////////////
    // Filter code ends here.
    ///////////////////////////////////////////

    // The final posterior is described by a cloud of weighted points.
    std::cout << "Fill histogram" << std::endl;
    for (ExampleSIR::SampleVector::iterator s = samples.begin();
             s != samples.end(); ++s) {
        hist1->Fill(s->second[0], s->second[1], s->first);
    }

    hist1->Draw("colz");
    gPad->Print("RunSimpleSIRExample.png");
}
#endif
