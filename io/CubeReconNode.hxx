#ifndef ReconNode_hxx_seen
#define ReconNode_hxx_seen

#include "TObject.h"

#include "CubeHandle.hxx"
#include "CubeReconObject.hxx"
#include "CubeReconState.hxx"
#include "CubeLog.hxx"

namespace Cube {
    class ReconNode;
    class ReconNodeContainer;
}

/// Class to contain the association between a ReconState and the object
/// assocated with the state.  This class is used to encode information
/// required to refit an object (e.g. refit a track as part of a global fit),
/// as well as to encode the path of extended objects.  This class can be
/// thought of as representing one point in a fit.
class Cube::ReconNode: public TObject {
public:
    ReconNode();
    virtual ~ReconNode();

    /// Get the state associated with this node.
    Cube::Handle<Cube::ReconState> GetState() const {return fState;}

    /// Set the state associated with this node.
    void SetState(Cube::Handle<Cube::ReconState>& state) {fState = state;}

    /// Get the reconstruction object associated with this node.
    Cube::Handle<Cube::ReconObject> GetObject() const {return fObject;}

    /// Set the reconstruction object associates with this node.
    void SetObject(Cube::Handle<Cube::ReconObject>& object) {fObject = object;}

    /// Get the goodness associated with the connection between the state and
    /// the object.
    double GetQuality() const {return fQuality;}

    /// Set the goodness associated with the connection between this state
    /// and object.
    void SetQuality(double quality) {fQuality = quality;}

    /// Print the object information.
    virtual void ls(Option_t *opt = "") const;

private:

    /// The state associated with the object.
    Cube::Handle<Cube::ReconState> fState;

    /// The object that is associated with the state.
    Cube::Handle<Cube::ReconObject> fObject;

    /// A log likelihood for the association of the object with the state.
    float fQuality;

    ClassDef(ReconNode,2);
};

/// A base class for containers of ReconNode objects.
class Cube::ReconNodeContainer
    : public TObject, public std::vector< Handle<ReconNode> > {
public:
    ReconNodeContainer();
    virtual ~ReconNodeContainer();

    /// Add a node to the container.  The std::vector::push_back method is
    /// overloaded so that it can check that the correct type of ReconState
    /// object is being added to the container.  This is overloaded by the
    /// ReconNodeContainerImpl template.
    virtual void push_back(const Cube::Handle<Cube::ReconNode>&) {
        throw std::runtime_error("Node container push_back directly called");
    }

    /// Print the object information.
    virtual void ls(Option_t *opt = "") const;


    ClassDef(ReconNodeContainer,1);
};

namespace Cube {
    /// Provide a class specific container for ReconNode objects.  This
    /// checks that all of the nodes contain the right class of state before
    /// they are added to the container.
    template <class T>
    class ReconNodeContainerImpl : public ReconNodeContainer {
    public:
        ReconNodeContainerImpl() {}
        virtual ~ReconNodeContainerImpl() {}

        /// Add a node to the container.  The std::vector::push_back method is
        /// overloaded so that it can make sure that the node contains the
        /// correct type of ReconState object.
        virtual void push_back(const Cube::Handle<Cube::ReconNode>& node) {

            // Check that node has a state
            if (node->GetState()){
                Cube::Handle<T> n(node->GetState());
                // Check that the state handle is valid (correct state type)
                if (!n) {
                    throw std::runtime_error(
                        "Wrong type of state being added to a "
                        "ReconNodeContainer");
                }
            }

            Cube::Handle<Cube::ReconObject> obj = node->GetObject();
            if (!obj) {
                throw std::runtime_error(
                    "Node added to a ReconNodeContainer"
                    " without an associated recon object");
            }
            std::vector< Cube::Handle<Cube::ReconNode> >::push_back(node);
        }

        ClassDefT(ReconNodeContainerImpl,1);
    };
    ClassDefT2(ReconNodeContainerImpl,T);
}
#endif
