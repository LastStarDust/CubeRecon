#include "CubeClusterManagement.hxx"

/// Check if two clusters meet the SFG definition of neighbors.
bool Cube::AreNeighbors(const Cube::ReconCluster& A,
                        const Cube::ReconCluster& B) {
    if (!(A.GetHitSelection())) return false;
    if (!(B.GetHitSelection())) return false;
    Cube::HitSelection::iterator b1 = A.GetHitSelection()->begin();
    Cube::HitSelection::iterator e1 = A.GetHitSelection()->end();
    Cube::HitSelection::iterator b2 = B.GetHitSelection()->begin();
    Cube::HitSelection::iterator e2 = B.GetHitSelection()->end();
    if (b1 == e1) return false;
    if (b2 == e2) return false;
    if ((*b1) == *(b2)) return true;
    if ((*b1) == (*(e2-1))) return true;
    if ((*(e1-1)) == *(b2)) return true;
    if ((*(e1-1)) == (*(e2-1))) return true;
    return false;
}
