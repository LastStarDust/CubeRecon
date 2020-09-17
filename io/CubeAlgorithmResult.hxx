#ifndef CubeAlgorithmResult_hxx_seen
#define CubeAlgorithmResult_hxx_seen

#include "CubeHandle.hxx"
#include "CubeReconObject.hxx"

#include <TNamed.h>

namespace Cube {
    class AlgorithmResult;
    class Algorithm;
}

/// A class to save results from a Algorithm object (this is the base class
/// to save the reconstruction results).  This contains information about the
/// Algorithm object that generated the results, as well as any user defined
/// results from running the algorithm.  The status of the algorithm is
/// returned by the AlgorithmResult::GetStatus() method.  The default status
/// is empty.  One significant use of the status field is to specify if an
/// algorithm has flagged an event to be saved or rejected.  If the status
/// contains "rej", then the event should be rejected.  The status can also be
/// used to specify which cuts an event has passed.  Typically, an algorithm
/// should create a status string, and added information using the "+="
/// operator.
///
/// By convention, the algorithm should fill the AlgorithmResult with some
/// standard Cube::ReconObjectContainers and Cube::HitSelections.  See the
/// AlgorithmResult documentation for more information.
///
/// * final -- A Cube::ReconObjectContainer that contains the main output
///            reconstruction objects for this algorithm.  Any Algorithms that
///            take this result as input will use the "final" reconstruction
///            objects.  This should be the last Cube:::ReconObjectContainer
///            added to the AlgorithmResult.
///
/// * used -- A Cube::HitSelection of the input hits that were used in this
///           event.  If the Algorithm did not create new hits, then the used
///           hit selection contains the final hits.  This should usually be
///           the last Cube::HitSelection added to the AlgorithmResult.  If
///           "used" is not the final hit selection, the correct hit selection
///           should be clearly documented.
///
/// * unused -- A Cube::HitSelection of input hits that were not used in this
///             event.
///
/// Other object containers and hit selections can be added to the
/// AlgorithmResult.  AlgorithmResult objects are often created by one
/// Algorithm, and then used as input to the next Algorithm.  In that case,
/// the Algorithm will use the last container of objects, and hit selection
/// for it's input.
class Cube::AlgorithmResult : public TNamed {
public:
    typedef std::vector<Cube::Handle<Cube::AlgorithmResult>> AlgorithmResults;
    typedef std::vector<Cube::Handle<Cube::ReconObjectContainer>> ReconObjects;
    typedef std::vector<Cube::Handle<Cube::HitSelection>> HitSelections;

    AlgorithmResult();
    explicit AlgorithmResult(const Algorithm& algo);
    AlgorithmResult(const char* name, const char* title);
    virtual ~AlgorithmResult();

    /// A static empty value for use as a default parameter in Algorithm.
    static const AlgorithmResult Empty;

    /// A copy constructor to turn a Cube::HitSelection object into a
    /// AlgorithmResult.  Through the magic of C++, this acts as a conversion
    /// operator from a Cube::HitSelection object and AlgorithmResult, so any
    /// routine that would accept a AlgorithmResult will also accept a
    /// Cube::HitSelection.
    AlgorithmResult(const Cube::HitSelection& hits);

    //@{
    /// Add a status string to the existing status.  The status string should
    /// be a single string without spaces, and it will be saved in the status
    /// string using a canonical format.  The format is best demonstrated by
    /// example:
    ///
    /// \code
    /// result->AddStatus("tpc");         // Save a flag
    /// result->AddStatus("good-track");  // Save another flag
    /// result->AddStatus("keep:true");   // Save a boolean value
    /// std::cout << result->GetStatus() << std::endl;
    /// \endcode
    ///
    /// Will print "(tpc) (good-track) (keep:true)".  Fields that are saved in
    /// the status can be verified with the following code idiom:
    ///
    /// \code
    /// if (result->GetStatus().find("(tpc)") != std::string::npos) {
    ///     // The TPC status field was found so do something...
    /// }
    /// \endcode
    ///
    void AddStatus(const char* status);
    void AddStatus(const std::string& status);
    //@}

    //@{
    /// Set a status summary of the algorithm result by over-writing the
    /// existing status.
    void SetStatus(const char* status);
    void SetStatus(const std::string& status);
    //@}

    /// Get the status summary of the result.  This can be used to summarize
    /// any error conditions that the algorithm found, cuts that were applied,
    /// or any other general status information.  The status summary is
    /// returned as a std::string and the default value is empty.  One
    /// typical use of this is to determine if an event should be saved or
    /// rejected by an analysis.
    std::string GetStatus(void) const;

    /// Add a new reconstruction object container to the AlgorithmResult.
    /// Also, the container will be owned by the AlgorithmResult.  The most
    /// recently added Cube::ReconObjectContainer object will become the default
    /// set of results.
    void AddObjectContainer(Cube::Handle<Cube::ReconObjectContainer> objects);

    /// Get a reconstruction object container out of this AlgorithmResult.
    /// The most recently added object container that matchs the requested
    /// name will be returned. This can be used to look for objects in
    /// sub-results by treating the name like a file system hierarchy
    /// (e.g. "algo1/subalgo1/subalgo2/objectName").
    Cube::Handle<Cube::ReconObjectContainer>
    GetObjectContainer(const char* object=NULL) const;

    /// Get the container of ReconObjectContainer objects. This is exposing
    /// some class internals, so you should usually prefer GetObjectContainer().
    const ReconObjects& GetObjectContainers() const {
        return fObjectContainers;}
    ReconObjects& GetObjectContainers() {return fObjectContainers;}

    /// Add a new hit selection into the AlgorithmResult.  Note that the hit
    /// selection name must be unique in this result.  Also, the hit selection
    /// will be owned by the AlgorithmResult.  The most recently added
    /// Cube::HitSelection object will become the default set of hits.
    void AddHitSelection(Cube::Handle<Cube::HitSelection> hits);

    /// Get a hit selection out of this AlgorithmResult.  The most recently
    /// added hit selection that matchs the requested name will be returned.
    Cube::Handle<Cube::HitSelection>
    GetHitSelection(const char* hit=NULL) const;

    /// Get the container of HitSelection objects.  This is exposing some
    /// class internals, so you should usually prefer GetHitSelection().
    const HitSelections& GetHitSelections() const {return fHitSelections;}
    HitSelections& GetHitSelections() {return fHitSelections;}

    /// Add a new sub AlgorithmResult.
    void AddAlgorithmResult(Cube::Handle<Cube::AlgorithmResult> result);

    /// Get a sub AlgorithmResult.  The most recently added algorithm result
    /// that matchs the requested name will be returned.  This can be used to
    /// look for sub-results by treating the name like a file system
    /// hierarchy (e.g. "algo1/subalgo1/subalgo2").
    Cube::Handle<Cube::AlgorithmResult>
    GetAlgorithmResult(const char* result=NULL) const;

    /// Get the container of children.  This is exposing some class internals,
    /// so you should prefer GetAlgorithmResult().
    const AlgorithmResults& GetResultsContainer() const {
        return fResultsContainer;
    }
    AlgorithmResults& GetResultsContainer() {return fResultsContainer;}

    /// Print the result information.
    virtual void ls(Option_t *opt = "") const;

    /// Return true if the algorithm is empty.
    bool IsEmpty() const;

protected:
    /// Initialize the object.
    void Initialize() {
        fStatusSummary.clear();
        fResultsContainer.clear();
        fObjectContainers.clear();
        fHitSelections.clear();
    }

    /// A summary of the algorithm status.  This is typical used to save the
    /// "event disposition" about whether an event is selected by a set of
    /// cuts or not, but can also be used to flag error conditions from the
    /// algorithm
    std::string fStatusSummary;

    /// A container of sub-algorithm results
    AlgorithmResults fResultsContainer;

    /// A container of reconstruction object containers.
    ReconObjects fObjectContainers;

    /// A container of hit selections.
    HitSelections fHitSelections;

    ClassDef(AlgorithmResult,3);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
