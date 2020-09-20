#include "CubeTreeRecon.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterHits.hxx"
#include "CubeSpanningTree.hxx"
#include "CubeFindKinks.hxx"
#include "CubeGrowClusters.hxx"
#include "CubeGrowTracks.hxx"
#include "CubeMergeXTalk.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>

Cube::TreeRecon::TreeRecon(): Cube::Algorithm("TreeRecon") { }

Cube::Handle<Cube::AlgorithmResult>
Cube::TreeRecon::Process(const Cube::AlgorithmResult& input,
                        const Cube::AlgorithmResult&,
                        const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Process TreeRecon" << std::endl;
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        CUBE_ERROR << "No hits" << std::endl;
        return result;
    }
    if (inputHits->empty()) {
        CUBE_ERROR << "Hits empty" << std::endl;
        return result;
    }

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // Apply all of the algorithms to the input.  The "do/while" construct
    // implements finalization.  The final objects from currentResult will be
    // copied into the final recon container.
    Cube::Handle<Cube::AlgorithmResult> currentResult;;
    do {
        // Cluster the 3D hits so they are all simply connected.
        Cube::Handle<Cube::AlgorithmResult> clusterHits
            = Run<Cube::ClusterHits>(*inputHits);
        if (!clusterHits) break;
        currentResult = clusterHits;
        result->AddAlgorithmResult(currentResult);

        // Further break the clusters based on the spanning tree.
        Cube::Handle<Cube::AlgorithmResult> spanningTree
            = Run<Cube::SpanningTree>(*currentResult);
        if (!spanningTree) break;
        currentResult = spanningTree;
        result->AddAlgorithmResult(currentResult);

        // Further break the clusters based on kinks.
        Cube::Handle<Cube::AlgorithmResult> findKinks
            = Run<Cube::FindKinks>(*currentResult);
        if (!findKinks) break;
        currentResult = findKinks;
        result->AddAlgorithmResult(currentResult);

        // Grow the clusters into track like objects.
        Cube::Handle<Cube::AlgorithmResult> growClusters
            = Run<Cube::GrowClusters>(*currentResult);
        if (!growClusters) break;
        currentResult = growClusters;
        result->AddAlgorithmResult(currentResult);

        // Grow the tracks to prevent gaps.
        Cube::Handle<Cube::AlgorithmResult> growTracks
            = Run<Cube::GrowTracks>(*currentResult);
        if (!growTracks) break;
        currentResult = growTracks;
        result->AddAlgorithmResult(currentResult);

        // Grow the tracks to prevent gaps.
        Cube::Handle<Cube::AlgorithmResult> mergeXTalk
            = Run<Cube::MergeXTalk>(*currentResult);
        if (!mergeXTalk) break;
        currentResult = mergeXTalk;
        result->AddAlgorithmResult(currentResult);

    } while (false);            // Always stops...

    // Copy the best result to the output
    if (currentResult) {
        Cube::Handle<Cube::ReconObjectContainer> objects
            = currentResult->GetObjectContainer("final");
        if (objects) {
            std::copy(objects->begin(), objects->end(),
                      std::back_inserter(*finalObjects));
        }
    }

    // Save the final objects last.
    result->AddObjectContainer(finalObjects);

    // Build the selections of used and unused hits.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
