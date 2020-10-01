#ifndef CubeInfo_hxx_seen
#define CubeInfo_hxx_seen

namespace Cube {
    class Info;
}

/// A collection of static functions and singleton classes that collect
/// detector related information.
class Cube::Info {
public:
    // Return the sub-detector for the sensor id.
    static int SubDetector(int i);

    /// @{ Return the number of the plane, bar or cube based on the sensor id.
    /// If the id is not for the cube detector, this will throw an exception.
    /// If the id is for a sensor, then the two relevant axis will be
    /// positive, and the sensor axis is -1.  If this is expanded for double
    /// ended readout, the "far" sensor will be "-2".  The "number" increments
    /// along the "X" axis.  The "bar" increments along the "Y" axis.  The
    /// "plane" increments along the "Z" axis.
    static int CubeNumber(int id);
    static int CubeBar(int id);
    static int CubePlane(int id);
    /// @}

    /// Return the projection for a hit based on the identifier.  The
    /// projections are X == 1 (0b001), Y == 2 (0b010), Z = 4 (0b100),
    /// YZ == 6 (0b110), XZ == 5 (0b101), XY == 3 (0b011),
    //  or XYZ == 7 (0b111).
    enum { kXAxis = 1, kYAxis = 2, kZAxis = 4,
           kYZProj = 6, kXZProj = 5, kXYProj = 3,
           kXYZProj = 7};

    static int IdentifierProjection(int id1);

    /// Build a cube identifer.
    static int Identifier3DST(int num, int bar, int pln);

    /// Take two or three hit identifiers and return the identifier of where
    /// the intersect.  If the ids do not correspond to a valid intersection,
    /// then return 0.
    static int Combine3DST(int id1, int id2, int id3 = -1);


};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
