#include <limits>
#include <cmath>
#include <algorithm>

#include "CubeCorrValues.hxx"
#include "CubeLog.hxx"

ClassImp(Cube::CorrValues);

/// Set the variance value for a free parameter.  The value is about 1E+154
/// (1E+19) for a double (float).
const Cube::CorrValues::Element Cube::CorrValues::kFreeValue =
    std::sqrt(std::numeric_limits<float>::max());
const Cube::CorrValues::Element Cube::CorrValues::kFreeThreshold =
    0.1*Cube::CorrValues::kFreeValue;
const Cube::CorrValues::Element
    Cube::CorrValues::kFixedValue = 1.0/Cube::CorrValues::kFreeValue;
const Cube::CorrValues::Element
    Cube::CorrValues::kFixedThreshold = 1.0/Cube::CorrValues::kFreeThreshold;

Cube::CorrValues::CorrValues(int n)
    : fVector(n), fMatrix(n), fNDOF(0), fTypeHash(0), fHessian(NULL) {
    // All parameters are initially free.
    for (int i=0; i<GetDimensions(); ++i) {
        fMatrix(i,i) = kFreeValue;
    }
}

Cube::CorrValues::CorrValues(double value, double uncertainty)
    : fVector(1), fMatrix(1), fNDOF(1), fTypeHash(0), fHessian(NULL) {
    SetValue(0,value);
    SetCovarianceValue(0,0,uncertainty*uncertainty);
}

Cube::CorrValues::CorrValues(const CorrValues& hv)
    : TObject(hv), fVector(hv.fVector), fMatrix(hv.fMatrix),
      fNDOF(hv.fNDOF), fTypeHash(hv.fTypeHash), fHessian(NULL) {}

Cube::CorrValues::CorrValues(const TVectorT<float>& v)
    : fVector(v), fMatrix(v.GetNoElements()),
      fNDOF(0), fTypeHash(0), fHessian(NULL) {
    // All parameters are initially free.
    for (int i=0; i<GetDimensions(); ++i) {
        fMatrix(i,i) = kFreeValue;
    }
}

Cube::CorrValues::CorrValues(const TVectorT<float>& v,
                               const TVectorT<float>& err)
    : fVector(v), fMatrix(v.GetNoElements()),
      fNDOF(0), fTypeHash(0), fHessian(NULL) {
    if (err.GetNoElements() != v.GetNoElements()) {
        CUBE_ERROR << "Mismatched number of elements" << std::endl;
        throw std::runtime_error("Invalid CorrValues constructor");
    }
    for (int i = 0; i < GetDimensions(); ++i) {
        fMatrix(i,i) = err(i);
    }
}

Cube::CorrValues::CorrValues(const TVectorT<float>& v,
                             const TMatrixTSym<float>& cov)
    : fVector(v), fMatrix(cov),
      fNDOF(0), fTypeHash(0), fHessian(NULL) {
    if (cov.GetNcols() != fVector.GetNoElements()) {
        CUBE_ERROR <<"Mismatch between element count and covariance columns."<<std::endl;
        throw std::runtime_error("Invalid CorrValues Constructor");
    }
    if (cov.GetNrows() != fVector.GetNoElements()) {
        CUBE_ERROR <<"Mismatch between element count and covariance rows."<<std::endl;
        throw std::runtime_error("Invalid CorrValues Constructor");
    }
}

void Cube::CorrValues::SetValues(const TVectorT<float>& v) {
    fVector = v;
}

void Cube::CorrValues::SetCovariance(const TMatrixTSym<float>& m) {
    fMatrix = m;
    fNDOF = 0;
    if (fHessian) delete fHessian;
    fHessian = NULL;
}

void Cube::CorrValues::SetValue(int i, double v) {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    fVector[i] = v;
}

double Cube::CorrValues::GetValue(int i) const {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    return fVector[i];
}

void Cube::CorrValues::SetCovarianceValue(int i, int j, double v) {
    if (i<0 || j<0) {
        CUBE_ERROR <<"Negative covariance index: (" << i << "," << j << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i || GetDimensions()<=j) {
        CUBE_ERROR <<"Negative covariance index: (" << i << "," << j << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    fMatrix(i,j) = v;
    fMatrix(j,i) = v;
    fNDOF = 0;
    if (fHessian) delete fHessian;
    fHessian = NULL;
}

double Cube::CorrValues::GetCovarianceValue(int i, int j) const {
    if (i<0 || j<0) {
        CUBE_ERROR <<"Negative covariance index: (" << i << "," << j << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i || GetDimensions()<=j) {
        CUBE_ERROR <<"Negative covariance index: (" << i << "," << j << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    return fMatrix(i,j);
}

double Cube::CorrValues::GetUncertainty(int i) const {
    double val = GetCovarianceValue(i,i);
    if (val>0) val = std::sqrt(val);
    return val;
}

void Cube::CorrValues::SetValue(double v) {
    if (GetDimensions() != 1) {
        CUBE_ERROR <<"Not 1 dimensional"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    SetValue(0,v);
}

double Cube::CorrValues::GetValue() const {
    if (GetDimensions() != 1) {
        CUBE_ERROR <<"Not 1 dimensional"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    return GetValue(0);
}

void Cube::CorrValues::SetUncertainty(double v) {
    if (GetDimensions() != 1) {
        CUBE_ERROR <<"Not 1 dimensional"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    SetCovarianceValue(0,0,v*v);
}

double Cube::CorrValues::GetUncertainty() const {
    if (GetDimensions() != 1) {
        CUBE_ERROR <<"Not 1 dimensional"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    return GetUncertainty(0);
}

void Cube::CorrValues::ResizeTo(int n) {
    fVector.ResizeTo(n);
    fMatrix.ResizeTo(n,n);
    fNDOF = 0;
    fTypeHash = 0;
    for (int i=0; i<GetDimensions(); ++i) {
        fMatrix(i,i) = kFreeValue;
    }
    if (fHessian) delete fHessian;
    fHessian = NULL;
}

int Cube::CorrValues::GetDimensions(void) const {
    return fVector.GetNoElements();
}

int Cube::CorrValues::GetNDOF() {
    if (fNDOF < 1) {
        fNDOF = 0;
        for (int i=0; i<GetDimensions(); ++i) {
            if (IsFixed(i)) continue;
            if (IsFree(i)) continue;
            ++ fNDOF;
        }
    }
    return fNDOF;
}

void Cube::CorrValues::SetType(const char* type) {
    // Use a modified Fowler, Noll, and Vo hash type 1a (FNV-1).  This is the
    // hash function uses the prime and offset suggested in TR1 for strings.
    // Note: The offset value for the normal FNV-1a is chosen so that an empty
    // string doesn't map to zero.  In TR1, the offset is chosen so that an
    // empty string will satisfy (hash_value == 0).
#define USE_TR1_FNV
#ifndef USE_TR1_FNV
    const int prime = 0x83; // The TR1 value
    const int offset = 0;   // The TR1 value
#else
    const int prime = 0x01000193;  // The FNV-1 value;
    const int offset = 0x811c9dc5; // The FNV-1 value;
#endif
    fTypeHash = offset;
    for (const char* c = type; *c != 0; ++c) {
        fTypeHash *= prime;
        fTypeHash ^= (unsigned int) *c;
    }
}

void Cube::CorrValues::SetFixed(int i) {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    for (int j=0; j<GetDimensions(); ++j) {
        fMatrix(j,i) = fMatrix(i,j) = 0.0;
    }
    fMatrix(i,i) = kFixedValue;
}

bool Cube::CorrValues::IsFixed(int i) const {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (fMatrix(i,i)<kFixedThreshold) return true;
    return false;
}

void Cube::CorrValues::SetFree(int i) {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    for (int j=0; j<GetDimensions(); ++j) {
        fMatrix(j,i) = fMatrix(i,j) = 0.0;
    }
    fMatrix(i,i) = kFreeValue;
}

bool Cube::CorrValues::IsFree(int i) const {
    if (i<0) {
        CUBE_ERROR <<"Negative element index: " << i<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    if (GetDimensions()<=i) {
        CUBE_ERROR <<"Out of bounds element index: " << i
                   << " (dim is " << GetDimensions() << ")"<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    return IsFree(fMatrix(i,i));
}

Cube::CorrValues Cube::CorrValues::Sum() const {
    double val = 0;
    double var = 0;

    for (int i=0; i<GetDimensions(); ++i) {
        val += GetValue(i);
        var += GetCovarianceValue(i,i);
        for (int j=i+1; j<GetDimensions(); ++j) {
            var += 2*GetCovarianceValue(i,j);
        }
    }

    return CorrValues(val,var);
}

void Cube::CorrValues::ls(Option_t *opt) const {
    fVector.ls(opt);
    fMatrix.ls(opt);
}

void Cube::CorrValues::Print(Option_t *opt) const {
    fVector.Print(opt);
    fMatrix.Print(opt);
}
bool Cube::CorrValues::Validate(bool fix) {
    int dim = GetDimensions();
    if (dim<1) return false;
    bool ok = true;
    for (int i=0; i < dim; ++i) {
        if (IsFixed(i)) SetFixed(i);
        if (IsFree(i)) SetFree(i);
    }
    for (int i=0; i < dim; ++i) {
        for (int j=i+1; j < dim; ++j) {
            // Check the off-diagonal elements to make sure they are equal.
            // I'm being overly careful here and checking both elements
            // against the average.
            Element mag = 0.5*(fMatrix(i,j) + fMatrix(j,i));
            if (std::abs(fMatrix(i,j)-mag)>2*fMatrix.GetTol()*mag) ok = false;
            if (std::abs(fMatrix(j,i)-mag)>2*fMatrix.GetTol()*mag) ok = false;
            if (fix) {
                fMatrix(i,j) = mag;
                fMatrix(j,i) = mag;
            }
        }
    }
    if (!ok && fix) {
        fNDOF = 0;
        if (fHessian) delete fHessian;
        fHessian = 0;
    }
    return ok;
}

const TMatrixTSym<float>& Cube::CorrValues::GetHessian() const {
    if (fHessian) return (*fHessian);
    fHessian = new TMatrixTSym<float>(GetDimensions());

    // Collect the non-zero elements of the covariance.
    std::vector<int> goodList;
    for (int i = 0; i<GetDimensions(); ++i) {
        if (IsFixed(i)) continue;
        if (IsFree(i)) continue;
        goodList.push_back(i);
    }

    if (goodList.size()>0) {
        // Make a symetric matrix of just the good elements and invert.
        TMatrixT<Element> h(goodList.size(),goodList.size());
        for (unsigned int i = 0; i<goodList.size(); ++i) {
            for (unsigned int j = i; j<goodList.size(); ++j) {
                h(j,i) = h(i,j) = fMatrix(goodList[i],goodList[j]);
            }
        }
        h.Invert();
        // Copy the values of the inverted covariance matrix into the Hessian.
        for (unsigned int i = 0; i<goodList.size(); ++i) {
            for (unsigned int j = 0; j<goodList.size(); ++j) {
                (*fHessian)(goodList[i],goodList[j]) = h(i,j);
            }
        }
    }

    // Set the Hessian values for the free and fixed variances.
    for (int i = 0; i<GetDimensions(); ++i) {
        if (IsFree(i)) (*fHessian)(i,i) = 1.0/kFreeValue;
        if (IsFixed(i)) (*fHessian)(i,i) = 1.0/kFixedValue;
    }
    return (*fHessian);
}

const Cube::CorrValues& Cube::CorrValues::operator =(const Cube::CorrValues&
                                                     rhs) {
    int dim = rhs.GetDimensions();

    if (dim != GetDimensions()) {
        ResizeTo(dim);
        if (fHessian) delete fHessian;
    }

    for (int i=0; i<dim; ++i) {
        fVector[i] = rhs.fVector[i];
        for (int j=0; j<dim; ++j) {
            fMatrix(i,j) = rhs.fMatrix(i,j);
        }
    }

    fNDOF = rhs.fNDOF;
    fTypeHash = rhs.fTypeHash;
    if (rhs.fHessian) fHessian = new TMatrixTSym<float>(*rhs.fHessian);

    return rhs;
}

Cube::CorrValues operator +(const Cube::CorrValues& a,
                           const Cube::CorrValues& b) {
    if (a.GetDimensions() != b.GetDimensions()) {
        CUBE_ERROR <<"Dimensions mismatch: "
                   << a.GetDimensions() << " != " << b.GetDimensions()<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    Cube::CorrValues result(a);

    int dim = a.GetDimensions();
    for (int i=0; i<dim; ++i) {
        result.SetValue(i,a.GetValue(i) + b.GetValue(i));
        for (int j=0; j<dim; ++j) {
            double val = a.GetCovarianceValue(i,j) + b.GetCovarianceValue(i,j);
            result.SetCovarianceValue(i,j,val);
        }
    }

    return result;
}

Cube::CorrValues operator +(double a, const Cube::CorrValues& b) {
    Cube::CorrValues result(b);
    int dim = b.GetDimensions();
    for (int i = 0; i<dim; ++i) result.SetValue(i,a + b.GetValue(i));
    return result;
}

Cube::CorrValues operator +(const Cube::CorrValues& a, double b) {return b + a;}

Cube::CorrValues operator -(const Cube::CorrValues& a,
                           const Cube::CorrValues& b) {
    if (a.GetDimensions() != b.GetDimensions()) {
        CUBE_ERROR <<"Dimensions mismatch: "
                   << a.GetDimensions() << " != " << b.GetDimensions()<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    Cube::CorrValues result(a);

    int dim = a.GetDimensions();
    for (int i=0; i<dim; ++i) {
        result.SetValue(i,a.GetValue(i) - b.GetValue(i));
        for (int j=0; j<dim; ++j) {
            double val = a.GetCovarianceValue(i,j) + b.GetCovarianceValue(i,j);
            result.SetCovarianceValue(i,j,val);
        }
    }

    return result;
}

Cube::CorrValues operator -(double a, const Cube::CorrValues& b) {
    Cube::CorrValues result(b);
    int dim = b.GetDimensions();
    for (int i = 0; i<dim; ++i) result.SetValue(i,a - b.GetValue(i));
    return result;
}

Cube::CorrValues operator -(const Cube::CorrValues& a, double b)  {
    Cube::CorrValues result(a);
    int dim = a.GetDimensions();
    for (int i = 0; i<dim; ++i) result.SetValue(i,a.GetValue(i)-b);
    return result;
}

Cube::CorrValues operator *(double a, const Cube::CorrValues& b) {
    Cube::CorrValues result(b);
    int dim = b.GetDimensions();
    for (int i=0; i<dim; ++i) {
        result.SetValue(i,a*b.GetValue(i));
        for (int j=0; j<dim; ++j) {
            result.SetCovarianceValue(i,j,a*b.GetCovarianceValue(i,j));
        }
    }
    return result;
}

Cube::CorrValues operator *(const Cube::CorrValues& b, double a) {
    return (a*b);
}

Cube::CorrValues operator /(const Cube::CorrValues& b, double a) {
    return ((1/a)*b);
}

Cube::CorrValues operator /(double a, const Cube::CorrValues& x) {
    Cube::CorrValues result(x);
    int dim = x.GetDimensions();
    for (int i=0; i<dim; ++i) {
        result.SetValue(i,a/x.GetValue(i));
        for (int j=0; j<dim; ++j) {
            double xi = x.GetValue(i);
            double xj = x.GetValue(j);
            double val = x.GetCovarianceValue(i,j)/(xi*xi*xj*xj);
            result.SetCovarianceValue(i,j,a*val);
        }
    }
    return result;
}

Cube::CorrValues operator *(const Cube::CorrValues& x,
                           const Cube::CorrValues& y) {
    if (x.GetDimensions() != y.GetDimensions()) {
        CUBE_ERROR <<"Dimensions mismatch: "
                   << x.GetDimensions() << " != " << y.GetDimensions()<<std::endl;
        throw std::runtime_error("CorrValues Range Error");
    }
    Cube::CorrValues result(x);

    int dim = y.GetDimensions();
    for (int i=0; i<dim; ++i) {
        result.SetValue(i,x.GetValue(i)*y.GetValue(i));
        for (int j=i; j<dim; ++j) {
            double cov = x.GetCovarianceValue(i,j)*y.GetValue(i)*y.GetValue(j);
            cov += y.GetCovarianceValue(i,j)*x.GetValue(i)*x.GetValue(j);
            if (cov<0) cov = -cov;
            result.SetCovarianceValue(i,j,cov);
            result.SetCovarianceValue(j,i,cov);
        }
    }
    return result;
}

Cube::CorrValues operator /(const Cube::CorrValues& x,
                           const Cube::CorrValues& y) {
    return x * (1.0/y);
}
