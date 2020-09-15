#ifndef CubeCompareReconObjects_hxx_seen
#define CubeCompareReconObjects_hxx_seen

#include <CubeReconCluster.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconVertex.hxx>
#include <CubeHandle.hxx>

namespace Cube {

    /// This is a predicate to order recon objects.  This puts the vertices,
    /// followed by tracks and clusters.  Within any type the objects are
    /// ordered by decreasing size (i.e. largest first).  If everything about
    /// two objects is equal, they are ordered by the value of the pointers.
    struct CompareReconObjects {
        bool operator () (Cube::Handle<Cube::ReconObject> lhs,
                          Cube::Handle<Cube::ReconObject> rhs) {
            Cube::Handle<Cube::ReconVertex> lv = lhs;
            Cube::Handle<Cube::ReconVertex> rv = rhs;
            Cube::Handle<Cube::ReconTrack> lt = lhs;
            Cube::Handle<Cube::ReconTrack> rt = rhs;
            Cube::Handle<Cube::ReconCluster> lc = lhs;
            Cube::Handle<Cube::ReconCluster> rc = rhs;

            if (lv && rv) return lv->GetNodes().size() > rv->GetNodes().size();
            if (lv && !rv) return true;
            if (!lv && rv) return false;

            if (lt && rt) return lt->GetNodes().size() > rt->GetNodes().size();
            if (lt && !rt) return true;
            if (!lt && rt) return false;

            if (lc && rc) {
                return (lc->GetHitSelection()->size()
                        > rc->GetHitSelection()->size());
            }
            if (lc && !rc) return true;
            if (!lc && rc) return false;

            return Cube::GetPointer(lhs) < Cube::GetPointer(rhs);
        }
    };
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
