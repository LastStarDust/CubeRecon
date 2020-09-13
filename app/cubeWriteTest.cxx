#include <iostream>

#include "CubeAlgorithmResult.hxx"

int main(int argc, char**argv) {
    std::cout << "cubeWriteTest" << std::endl;

    Cube::Handle<Cube::AlgorithmResult> inputResult(
        new Cube::AlgorithmResult("InputHits"));

    Cube::Handle<Cube::HitSelection> testHits(
        new Cube::HitSelection("testHits"));

    Cube::Handle<Cube::AlgorithmResult> testResult(
        new Cube::AlgorithmResult(*testHits));

    return 0;
}
