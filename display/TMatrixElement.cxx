#include "TMatrixElement.hxx"

#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoEltu.h>
#include <TGeoMatrix.h>
#include <TVectorD.h>

Cube::TMatrixElement::~TMatrixElement() {}

Cube::TMatrixElement::TMatrixElement(const char* name,
                                   const TVector3& position,
                                   const TMatrixD& matrix,
                                   bool longAxis)
    : TEveGeoShape(name) {
    Initialize(position,matrix,longAxis);
}

Cube::TMatrixElement::TMatrixElement(const char* name,
                                   const TVector3& position,
                                   const TMatrixF& matrix,
                                   bool longAxis)
    : TEveGeoShape(name) {
    TMatrixD mat(matrix);
    Initialize(position,mat,longAxis);
}

void Cube::TMatrixElement::Initialize(const TVector3& position,
                                      const TMatrixD& matrix,
                                      bool longAxis) {
    // Find the rotation of the object to be displayed.  If longAxis is true,
    // the the matrix is represented as a tube with the long axis along the
    // local Z direction, and the major and minor in the XY plane.  Otherwise,
    // the matrix is represented as a tube with the short axis along the local
    // Z direction.
    TMatrixD tubeRot(3,3);
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            tubeRot(i,j) = matrix(i,j);
            if (!std::isfinite(matrix(i,j))) {
                fValid = false;
                return;
            }
        }
    }

    double tubeAxis;
    double tubeMajor;
    double tubeMinor;
    TVectorD values(3);
    TMatrixD tubeEigen(tubeRot.EigenVectors(values));

    // Check that the eigenvalues are real.
    for (int i=0; i<3; ++i) {
        if (!std::isfinite(values(i))) {
            fValid = false;
            return;
        }
        for (int j=0; j<3; ++j) {
            if (!std::isfinite(tubeEigen(i,j))) {
                fValid = false;
                return;
            }
        }
    }

    if (longAxis) {
        tubeAxis = std::sqrt(values(0));
        tubeRot(0,2) = tubeEigen(0,0);
        tubeRot(1,2) = tubeEigen(1,0);
        tubeRot(2,2) = tubeEigen(2,0);
        tubeMajor = std::sqrt(values(1));
        tubeRot(0,0) = tubeEigen(0,1);
        tubeRot(1,0) = tubeEigen(1,1);
        tubeRot(2,0) = tubeEigen(2,1);
        tubeMinor = std::sqrt(values(2));
        tubeRot(0,1) = tubeEigen(0,2);
        tubeRot(1,1) = tubeEigen(1,2);
        tubeRot(2,1) = tubeEigen(2,2);
    }
    else {
        tubeAxis = std::sqrt(values(2));
        tubeRot(0,2) = tubeEigen(0,2);
        tubeRot(1,2) = tubeEigen(1,2);
        tubeRot(2,2) = tubeEigen(2,2);
        tubeMajor = std::sqrt(values(0));
        tubeRot(0,0) = tubeEigen(0,0);
        tubeRot(1,0) = tubeEigen(1,0);
        tubeRot(2,0) = tubeEigen(2,0);
        tubeMinor = std::sqrt(values(1));
        tubeRot(0,1) = tubeEigen(0,1);
        tubeRot(1,1) = tubeEigen(1,1);
        tubeRot(2,1) = tubeEigen(2,1);
    }

    // Create the rotation matrix.
    TGeoRotation rot;
    rot.SetMatrix(tubeRot.GetMatrixArray());

    // Set the translation
    TGeoTranslation trans(position.X(), position.Y(), position.Z());

    // Finally set the transform for the object.
    TGeoCombiTrans rotTrans(trans,rot);
    SetTransMatrix(rotTrans);

    // Make sure the tube size doesn't get too small.
    tubeAxis = std::max(1.5, tubeAxis);
    tubeMajor = std::max(1.5, tubeMajor);
    tubeMinor = std::max(1.5, tubeMinor);

    // Create the shape to display.  This has to play some fancy footsie to
    // get the gGeoManager memory management right.  It first saves the
    // current manager, then gets an internal geometry manager used by
    // TEveGeoShape, and then resets the old manager once the shape is
    // created.  You gotta love global variables...
    TGeoManager* saveGeom = gGeoManager;
    gGeoManager = GetGeoMangeur();
    TGeoShape* geoShape = new TGeoEltu(tubeMajor,tubeMinor,tubeAxis);
    SetShape(geoShape);
    gGeoManager = saveGeom;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
