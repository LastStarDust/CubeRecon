#ifndef CubeCorrValues_hxx_seen
#define CubeCorrValues_hxx_seen

#include <string>
#include <map>
#include <vector>
#include <typeinfo>
#include <algorithm>
#include <iostream>

#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>

namespace Cube {
    class CorrValues;
}

/// A class representing a vector of values and the correlations between them.
/// The CorrValues class is primarily used to represent experimental values
/// and their covariance matrix. The covariance matrix represents a free
/// parameter by setting a very large variance, and a fixed parameter by
/// setting a very small variance (approximately the maximum and minimum
/// floating point numbers).
///
/// The CorrValues class defines several operations to simplify the
/// propagation of errors, and is particularly optimized for the uncertainty
/// of one variable.  For example if you have two parameters with
/// uncertainty.
///
/// \code
/// X +- sigmaX
/// Y +- sigmaY
/// \endcode
///
/// where sigmaX and sigmaY is the standard uncertainty, then you can use
/// CorrValues objects to do arithmetic with X and Y.
///
/// \code
/// CorrValues Xvalue(X,sigmaX);
/// CorrValues Yvalue(Y,sigmaY);

/// CorrValues XYsum = Xvalue + Yvalue; // The sum of X and Y
/// CorrValues XYprod = Xvalue*Yvalue;  // The product of X and Y
/// CorrValues XYquot = Xvalue/Yvalue;  // The quotient of X and Y
/// \endcode
///
/// Operations with constant values (doubles) are also defined
///
/// \code
/// double A;
/// CorrValues Xvalue(X,sigmaX);
///
/// CorrValues AXsum = A + Xvalue;  // The sum of A and X
/// CorrValues AXprod = A*Xvalue;   // The product of A and X
/// CorrValues XAquot = Xvalue/A;   // The quotient of X and A
/// \endcode
///
/// These operations are defined for any dimension of CorrValues.  For
/// dimensions larger than one, the full covariance matrix is propagated.
class Cube::CorrValues: public TObject {
    /// The type of vector representation.  The code assumes that one of the
    /// TVectorT<> types is used.
    typedef TVectorT<double> Vector;

    /// The type of matrix representation.  The code assumes that one of the
    /// TMatrixT<> types is used.  The matrix should always be symmetric.
    typedef TMatrixTSym<double> Matrix;

public:
    /// The type of floating point number to use in a class.
    typedef float Element;

    /// Construct an empty n dimensional object [default: 0 dimension].
    CorrValues(int n=0);

    /// Construct a one dimensional parameter with an uncertainty.  The
    /// variance is (as expected) the uncertainty squared.
    CorrValues(double value, double uncertainty);

    /// The copy constructor.
    CorrValues(const CorrValues& hv);

    /// Construct an object from a vector of parameters, but without a
    /// matrix.  The matrix is taken to be zero meaning that all values are
    /// fixed.
    CorrValues(const TVectorT<float> &v);

    /// Construct an object from a vector of parameters, and a vector of
    /// uncertainties.  The uncertainties are the square-root of the
    /// variances, and are used to construct a diagonal matrix.
    CorrValues(const TVectorT<float> &v,const TVectorT<float> &vE);

    /// Construct an object from a vector of parameters and a matrix.  The
    /// parameter vector and matrix must have the same number of dimensions.
    CorrValues(const TVectorT<float> &v,const TMatrixTSym<float> &C);

    virtual ~CorrValues(){};

    /// The assignment operator
    const CorrValues& operator =(const CorrValues& rhs);

    /// Assign the values in the vector.
    void SetValues(const TVectorT<float>& v);

    /// Return the vector of parameters.
    const TVectorT<float>& GetValues() const {return fVector;}

    /// Assign the values in the matrix
    void SetCovariance(const TMatrixTSym<float>& m);

    /// Return the matrix of correlations.
    const TMatrixTSym<float>& GetCovariance() const {return fMatrix;}

    /// Set the value of a parameter.
    void SetValue(int i, double v);

    /// Get the value of a parameter;
    double GetValue(int i) const;

    /// Set the covariance of a parameter.
    void SetCovarianceValue(int i, int j, double v);

    /// Get the covariance of a parameter.
    double GetCovarianceValue(int i, int j) const;

    /// Get the uncertainty (standard error) on a parameter.  This does not
    /// take into account correlations with other parameters.
    double GetUncertainty(int i) const;

    /// @{Provide convenience methods to access the value, and uncertainty for
    /// one parameter CorrValues objects.
    void SetValue(double v);
    double GetValue() const;
    void SetUncertainty(double v);
    double GetUncertainty() const;
    /// @}

    /// Get the Hessian matrix (the inverse of the covariance), also known as
    /// the curvature matrix.  If the covariance as free or fixed parameters,
    /// the corresponding rows and columns in the hessian will be zero.
    const TMatrixTSym<float>& GetHessian() const;

    /// Set the dimensionality of the object.  This has the side effect
    /// of clearing the object.
    void ResizeTo(int n);

    /// Return the number of dimensions in the object.
    int GetDimensions() const;

    /// Return the number of degrees of freedom in the object.  This is
    /// usually the same as the number of object dimensions.
    int GetNDOF();

    /// Set the number of degrees of freedom.  This will usually be the same
    /// as the number of object dimensions, but may be reduced if some of
    /// the elements have an infinite uncertainty.  This doesn't affect the
    /// object calculations, but can be used to pass the information
    /// along.
    void SetNDOF(int ndf) {fNDOF = ndf;}

    /// Get the object type hash.  The type hash reflects the definitions
    /// of the object parameters.  It can be used to help keep
    /// incompatible objects from being used in the same expression.
    unsigned int GetTypeHash() const  {return fTypeHash;}

    /// Set the object type has.  The hash is constructed of the field
    /// definitions.
    void SetType(const char* type);

    /// Mark that parameter is fixed.
    void SetFixed(int param);

    /// Check to see if a parameter is fixed.
    bool IsFixed(int param) const;

    /// Check to see if a variance value means that the parameter is fixed.
    static bool IsFixed(double variance) {return (variance<kFixedThreshold);}

    /// Mark that a parameter is free.
    void SetFree(int param);

    /// Check to see if a parameter is free.
    bool IsFree(int param) const;

    /// Check to see if a variance value means that the parameter is free.
    static bool IsFree(double variance) {return (kFreeThreshold<variance);}

    /// The sum of the components as a one-dimensional correlated value.
    CorrValues Sum() const;

    /// List the correlated values information.
    virtual void ls(Option_t *opt = "") const;

    /// Print the contents of the correlated values.
    virtual void Print(Option_t *opt = "") const;

    /// Make sure that the covariance matrix is valid.  If the argument is
    /// true, then this ensures that the matrix is symmetric.
    bool Validate(bool fix=true);

    /// The value assigned to a variance for a free paramter.  This must be
    /// larger than kFreeThreshold.
    static const Element kFreeValue;

    /// The value assigned to a variance for a fixed paramter.  This must be
    /// smaller than kFixedThreshold.
    static const Element kFixedValue;

private:
    /// If a variance is larger than kFreeThreshold, the parameter is
    /// considered to be free.
    static const Element kFreeThreshold;

    /// If a variance is smaller than kFixedThreshold, the parameter is
    /// considered to be fixed.
    static const Element kFixedThreshold;

    /// The vector of parameter values.
    TVectorT<float> fVector;

    /// The vector of parameer correlations.
    TMatrixTSym<float> fMatrix;

    /// The number of independent parameters.
    int fNDOF;

    /// A semi-unique hash of the object element definitions.
    unsigned int fTypeHash;

    /// The value of the Hessian (if it's already calculated).
    mutable TMatrixTSym<float>* fHessian; //! Don't Save

    ClassDef(CorrValues,3);
};

///@{ Define basic additive arithmetic operations.  These operations are
///defined similar to vector addition so that for example,
/// \code
/// X+Y == [x1,x2]+[y1,y2] --> [x1+y1,x2+y2]
/// cov(X+Y) = cov(x1+y1, x2+y2) = cov(x1,x2) + cov(y1,y2)
/// \endcode
/// Addition by a constant A is defined as
/// \code
/// A+X = [x1, x2] + (a,a) --> [X1 + A, X2 + A]
/// cov(A+X) = cov(A+x1, A+x2) = cov(x1,x2)
/// \endcode
Cube::CorrValues operator +(const Cube::CorrValues& x, const Cube::CorrValues& y);
Cube::CorrValues operator +(double a, const Cube::CorrValues& x);
Cube::CorrValues operator +(const Cube::CorrValues& x, double a);
Cube::CorrValues operator -(const Cube::CorrValues& x, const Cube::CorrValues& y);
Cube::CorrValues operator -(double a, const Cube::CorrValues& x);
Cube::CorrValues operator -(const Cube::CorrValues& x, double a);
/// @}

/// Multiplication of CorrValues.  This is most useful for the multiplication
/// of one parameter correlated values, but is defined for any dimensionality.
/// For multi-dimensional correlated values, this is the "zip" of the vectors.
/// \code
/// X*Y == [x1,x2]*[y1,y2] --> [x1*y1,x2*y2]
/// \endcode
Cube::CorrValues operator *(const Cube::CorrValues& x, const Cube::CorrValues& y);

/// Division of a constant by a set of correlated values.  This is defined as
/// \code
/// A/X = A/[X1, X2] --> [A/X1, A/X2]
/// cov(A/X) = cov(A/X1, A/X2) = A*cov(1/X1, 1/X2)
/// \endcode
Cube::CorrValues operator /(double a, const Cube::CorrValues& x);

/// Division of CorrValues.  This is defined in terms of
/// \code
/// X/Y = X*(1.0/Y)
/// \endcode
/// and yields a correlated value of dimension one.
Cube::CorrValues operator /(const Cube::CorrValues& x, const Cube::CorrValues& y);

/// @{ Define basic multiplication by a constant.  These operators are defined
/// similar to vector multiplication by a constant.
/// \code
/// A*X = A*[x1, x2]  --> [A*X1, A*X2]
/// cov(A*X) = cov(A*x1, A*x2) = A*cov(x1,x2)
/// \endcode
Cube::CorrValues operator *(double a, const Cube::CorrValues& x);
Cube::CorrValues operator *(const Cube::CorrValues& x, double a);
Cube::CorrValues operator /(const Cube::CorrValues& x, double a);
/// @}
#endif
