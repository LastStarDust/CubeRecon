#ifndef CubeClusterManagement_hxx_seen
#define CubeClusterManagement_hxx_seen

#include <CubeReconCluster.hxx>

#include <algorithm>

namespace Cube {

    /// This takes iterators for two ordered collections of Cube hits and
    /// combines them into a new ordered collection.  The iterators are
    /// generally for hits (so Cube::HitSelection::iterator), but any
    /// iterator that provides a "pointer" to an object with a method
    /// "TVector GetPosition()" will work.  Assuming that clusters A and B
    /// are found to be neighbors (which means they share a
    /// beginning/ending hit), and that they should be combined, then
    ///
    /// \code
    /// Cube::HitSelection output;
    /// Cube::HitSelection split
    ///     = Cube::CombineNeighbors(
    ///           A->GetHitSelection()->begin(),A->GetHitSelection->end(),
    ///           B->GetHitSelection()->begin(),B->GetHitSelection->end(),
    ///           output);
    /// \endcode
    ///
    /// Will put the hits for A and B into the output hit selection in the
    /// right order.  The hits from the first cluster will always be
    /// first, but may be reversed from the original (the ending of the
    /// first cluster may be the beginning of the output cluster).  The
    /// hits from the second cluster are at the end of the output cluster.
    /// Because the clusters must share a single hit, the output cluster
    /// will have one less than the total number of hits in the two input
    /// clusters.  The return value is the iterator for the hit where the
    /// two clusters joined.  That hit will have been in both clusters.
    template<typename iterator, typename container>
    typename container::iterator CombineNeighbors(iterator b1, iterator e1,
                                                  iterator b2, iterator e2,
                                                  container& output);

    /// This returns true if two clusters share a hit at either end of the
    /// cluster.  This is treating a cluster as an ordered list of hits (which
    /// is a perversion, but works OK here).
    bool AreNeighbors(const Cube::ReconCluster& A,
                      const Cube::ReconCluster& B);


    /// Create clusters with one hit per cluster.  The clusters are ordered
    /// according to the input iterators.
    template<typename iterator>
    Cube::Handle<Cube::ReconObjectContainer>
    CreateHitClusters(iterator begin, iterator end);

    /// Create a cluster with the correct covariance from a set of. The
    /// iterators should resolve to Cube::Handle<Cube::Hit> object.
    template<typename iterator>
    Cube::Handle<Cube::ReconCluster>
    CreateCluster(const char* name, iterator begin, iterator end);
}

template<typename iterator, typename container>
typename container::iterator
Cube::CombineNeighbors(iterator b1, iterator e1,
                       iterator b2, iterator e2,
                       container& output) {
    typename container::iterator middle;
    if (b1 == b2 && e1 == e2) {
        CUBE_ERROR << "Cannot combine a cluster with a twin" << std::endl;
        return output.end();
    }
    int maxSize = (e1-b1) + (e2-b2);
    if (maxSize < 1) {
        CUBE_ERROR << "Inconceivable!" << std::endl;
        return output.end();
    }
    output.clear();
    output.reserve(maxSize);
    if ((*(e1-1)) == (*b2)) {
        std::copy(b1,e1,std::back_inserter(output));
        middle = output.end()-1;
        std::copy(b2+1,e2,std::back_inserter(output));
        return middle;
    }
    if ((*(e1-1)) == (*(e2-1))) {
        std::copy(b1,e1,std::back_inserter(output));
        middle = output.end()-1;
        --e2; while (e2 != b2) *(std::back_inserter(output)) = *(--e2);
        return middle;
    }
    if ((*b1) == *(b2)) {
        while (e1 != b1) *(std::back_inserter(output)) = *(--e1);
        middle = output.end()-1;
        std::copy(b2+1,e2,std::back_inserter(output));
        return middle;
    }
    if ((*b1) == (*(e2-1))) {
        while (e1 != b1) *(std::back_inserter(output)) = *(--e1);
        middle = output.end()-1;
        --e2; while (e2 != b2) *(std::back_inserter(output)) = *(--e2);
        return middle;
    }
    CUBE_ERROR << "Trying to combine clusters that are not neighbors"
               << std::endl;
    return output.end();
}

/// Take hits in a container between begin and end, and create one cluster per
/// hit.
template<typename iterator>
Cube::Handle<Cube::ReconObjectContainer>
Cube::CreateHitClusters(iterator begin, iterator end) {
    Cube::Handle<Cube::ReconObjectContainer> out(
        new Cube::ReconObjectContainer);
    while (begin != end) {
        Cube::Handle<Cube::ReconCluster> cluster
            = Cube::CreateCluster("cluster", begin,begin+1);
        out->push_back(cluster);
        ++begin;
    }
    return out;
}

/// Take some hits and create a cluster.  This uses
/// Cube::ReconCluster::FillFromHits.  The return value is a handle to the new
/// cluster.  The hits in the cluster are ordered according to the input
/// iterators.
template<typename iterator>
Cube::Handle<Cube::ReconCluster>
Cube::CreateCluster(const char* name, iterator begin, iterator end) {
    if (begin == end) return Cube::Handle<Cube::ReconCluster> ();
    Cube::Handle<Cube::ReconCluster> cluster(new Cube::ReconCluster);
    cluster->FillFromHits(name,begin,end);
    // Make an estimate of the effective number of hits in the cluster.
    // Large charge hits count more than small charge hits.
    double hits = 0.0;
    double avg = 0.0;
    // The scale, with a value of 4, gives an approximation of the relative
    // variance on the moments for each hit.  The effective number of hits in
    // the moment is calculated using this.
    double scale = 4.0;
    for (Cube::HitSelection::iterator h = cluster->GetHitSelection()->begin();
         h != cluster->GetHitSelection()->end(); ++h) {
        hits += (*h)->GetCharge();
        avg += std::pow((*h)->GetCharge(),scale);
    }
    avg /= cluster->GetHitSelection()->size();
    avg = std::pow(avg,1.0/scale);
    hits /= avg;
    // Fix the uncertainty of the cluster so that it's right for a single hit,
    // and reasonable for clusters including crosstalk.
    double xx = cluster->GetMoments()(0,0)/hits;
    double yy = cluster->GetMoments()(1,1)/hits;
    double zz = cluster->GetMoments()(2,2)/hits;
    double tt = cluster->GetPositionVariance().T();
    Cube::Handle<Cube::ClusterState> state = cluster->GetState();
    state->SetPositionVariance(xx,yy,zz,tt);
    return cluster;
}
#endif
