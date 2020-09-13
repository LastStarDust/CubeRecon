#ifndef CubeAlgorithm_hxx_seen
#define CubeAlgorithm_hxx_seen

#include <string>

#include "CubeHitSelection.hxx"
#include "CubeAlgorithmResult.hxx"

namespace Cube {
    class Algorithm;
}

/// Base class for reconstruction algorithms.  It contains the basic
/// information about the "flavor" and version.  Specific algorithms should
/// inherit the class and returned results using a AlgorithmResult object.
/// Inherited classes must implement the method
///
/// \code
/// Cube::Handle<AlgorithmResult>
/// Process(const Cube::AlgorithmResult& input
///         const Cube::AlgorithmResult& input1
///                       = Cube::AlgorithmResult::Empty,
///         const Cube::AlgorithmResult& input2
///                       = Cube::AlgorithmResult::Empty)
/// \endcode
///
/// which is a pure virtual member of this class.  The user should use
/// Algorithm::CreateResult() to construct an empty result that is returned
/// by the Process method.  This convenience function will set the algorithm
/// result name, and fill the bookkeeping information.
///
/// By convention, the algorithm should fill the AlgorithmResult with some
/// standard TReconObjectContainers and THitSelections.  See the
/// AlgorithmResult documentation for more information.
///
/// AlgorithmResult objects are often created by one Algorithm, and then
/// used as input to the next Algorithm.  In that case, the Algorithm will
/// use the last container of objects, and hit selection for it's input.
///
class Cube::Algorithm : public TNamed {
public:
    Algorithm(const char* name, const char* title="An Algorithm");
    virtual ~Algorithm();

    /// This method returns information about the compiled version
    virtual const std::string& GetVersion() const {return fVersion;}

    /// The routine that does the actual work.  This must be implemented by
    /// any derived class.  By convention, the derived class should only name
    /// the parameters that are actually used.  For instance, a derived class
    /// that only uses the first AlgorithmResult would be declared as
    /// \code
    /// Cube::Handle<Cube::AlgorithmResult>
    /// Process(const Cube::AlgorithmResult& input,
    ///         const Cube::AlgorithmResult& input1
    ///                    = Cube::AlgorithmResult::Empty,
    ///         const Cube::AlgorithmResult& input2
    ///                    = Cube::AlgorithmResult::Empty);
    /// \endcode
    /// and defined as
    /// \code
    /// Cube::Handle<Cube::AlgorithmResult>
    /// Process(const Cube::AlgorithmResult& input,
    ///         const Cube::AlgorithmResult&,
    ///         const Cube::AlgorithmResult&) { ... }
    /// \endcode
    virtual Cube::Handle<Cube::AlgorithmResult>
    Process(const Cube::AlgorithmResult& in0,
            const Cube::AlgorithmResult& in1 = Cube::AlgorithmResult::Empty,
            const Cube::AlgorithmResult& in2 = Cube::AlgorithmResult::Empty
        ) = 0;

    /// Reset algorithm containers, but keep do not change the
    /// parameters by default
    virtual void Clear(Option_t* option);

protected:
    /// Set the version of the algorithm.  This should be set in the derived
    /// classes constructor.
    void SetVersion(const char* v);

    /// A convenient way to create an empty algorithm result to be returned by
    /// a derived class.
    Cube::Handle<Cube::AlgorithmResult> CreateResult() {
        return Cube::Handle<AlgorithmResult>(new AlgorithmResult(*this));
    }

    /// Templates to simplify calling sub-algorithms.  These handle the
    /// Algorithm memory management.
    template<typename T>
    Cube::Handle<Cube::AlgorithmResult> Run(const Cube::AlgorithmResult& in1) {
        std::unique_ptr<T> ptr(new T);
        return ptr->Process(in1);
    }

    template<typename T>
    Cube::Handle<Cube::AlgorithmResult> Run(const Cube::AlgorithmResult& in1,
                                            const Cube::AlgorithmResult& in2) {
        std::unique_ptr<T> ptr(new T);
        return ptr->Process(in1,in2);
    }

    template<typename T>
    Cube::Handle<Cube::AlgorithmResult> Run(const Cube::AlgorithmResult& in1,
                                            const Cube::AlgorithmResult& in2,
                                            const Cube::AlgorithmResult& in3) {
        std::unique_ptr<T> ptr(new T);
        return ptr->Process(in1,in2,in3);
    }

private:
    /// A version string that must be set in the constructor with SetVersion().
    std::string fVersion;

    ClassDef(Algorithm, 1);
};
#endif
