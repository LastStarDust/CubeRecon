#ifndef CubeSafeLine_hxx_seen
#define CubeSafeLine_hxx_seen

#include <TVector3.h>
#include <TPrincipal.h>

#include <memory>

namespace Cube {
    template<typename iterator>
    bool SafeLine(iterator begin, iterator end,
                  TVector3& pnt, TVector3& dir);

    template<typename iterator>
    bool SafeLine(iterator begin, iterator end,
                  TVector3& pnt, TVector3& dir,
                  double& chi2, double& ndof);

    template<typename iterator>
    double SafeChi2(iterator begin, iterator end,
                    const TVector3& pnt, const TVector3& dir,
                    double& ndof);
}

/// Templates to take Cube hits and safely fit a line to them in 3D.  It's
/// playing tricks so that there is a reasonable answer as long as there are
/// at least two hits.  The iterators are generally for hits (so
/// ND::THitSelection::iterator), but any iterator that provides a "pointer"
/// to an object with a method "TVector GetPosition()" will work.  This
/// returns true as long as there isn't an error.
template<typename iterator>
bool Cube::SafeLine(iterator begin, iterator end,
                       TVector3& pnt, TVector3& dir) {
    if ((end-begin) < 2) {
        return false;
    }

    // Use principal component analysis to estimate the line.
    TPrincipal pca(3,"");
    for (iterator h = begin; h != end; ++h) {
        double row[3] =  {(*h)->GetPosition().X(),
                          (*h)->GetPosition().Y(),
                          (*h)->GetPosition().Z()};
        // Should this be charge weighting??? There are ways...
        pca.AddRow(row);
    }
    pca.MakePrincipals();

    // Find the point.
    double pcaPoint[3] = {0.0, 0.0, 0.0};
    double point[3];
    pca.P2X(pcaPoint,point,3);
    pnt = TVector3(point);

    pcaPoint[0] = 1.0;
    pca.P2X(pcaPoint,point,3);
    dir = TVector3(point) - pnt;

    if (dir.Mag() < 0.01) return false;
    dir = dir.Unit();

    return true;
}

/// A template to calculate the chi2 for hits to a line defined by a point and
/// a direction.  It returns the chi2, and the number of degrees of freedom.
/// If there aren't enough points (i.e. at least 3), it returns zero degrees
/// of freedom.  For two points, the chi2 might still be non zero.
template<typename iterator>
double Cube::SafeChi2(iterator begin, iterator end,
                         const TVector3& pnt, const TVector3& dir,
                         double& ndof) {
    if ((end-begin) < 2) {
        ndof = 0.0;
        return 0;
    }

    // Set the number of degrees of freedom.  There are six things determined
    // (xyz, dxdydx), but only 4 of them are independent.  The degrees of
    // freedom for the hit are the distances to the line (in 2D).
    ndof = 2.0*(end-begin) - 4.0;

    // Find the chi2.
    double variance = 100.0/12.0;
    double chi2 = 0.0;
    for (iterator h = begin; h != end; ++h) {
        TVector3 offset = (*h)->GetPosition() - pnt;
        double shift = offset * dir;
        offset = offset - shift*dir;
        chi2 += offset.Mag2()/variance;
    }

    return chi2;
}


/// A template to take Cube hits, safely fit a line to them in 3D, and
/// calculate the chi2/ndof.  It's playing tricks so that there is a
/// reasonable answer as long as there are at least two hits.  This returns
/// true as long as there isn't an error.
template<typename iterator>
bool Cube::SafeLine(iterator begin, iterator end,
                       TVector3& pnt, TVector3& dir,
                       double& chi2, double& ndof) {
    chi2 = 0.0;
    ndof = 0.0;
    if (!SafeLine(begin,end,pnt,dir)) return false;
    chi2 = SafeChi2(begin,end,pnt,dir,ndof);
    return true;
}
#endif
