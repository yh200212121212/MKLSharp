#pragma once

using namespace System;
using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Lapack {
  public:
    static __int64 sgetrf(LapackLayout Layout, int m, int n,
                          array<float>^ a, int lda, [Out]array<__int64>^% ipiv);
    static __int64 dgetrf(LapackLayout Layout, int m, int n,
                          array<double>^ a, int lda, [Out]array<__int64>^% ipiv);

    static __int64 sgetrs(LapackLayout Layout, char trans,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<__int64>^ ipiv, array<float>^ b, int ldb);
    static __int64 dgetrs(LapackLayout Layout, char trans,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<__int64>^ ipiv, array<double>^ b, int ldb);

    static __int64 sgeequ(LapackLayout Layout, int m, int n,
                          array<float>^ a, int lda,
                          [Out]array<float>^% r, [Out]array<float>^% c,
                          [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax);
    static __int64 dgeequ(LapackLayout Layout, int m, int n,
                          array<double>^ a, int lda,
                          [Out]array<double>^% r, [Out]array<double>^% c,
                          [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax);

    static __int64 sgeequb(LapackLayout Layout, int m, int n,
                           array<float>^ a, int lda,
                           [Out]array<float>^% r, [Out]array<float>^% c,
                           [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax);
    static __int64 dgeequb(LapackLayout Layout, int m, int n,
                           array<double>^ a, int lda,
                           [Out]array<double>^% r, [Out]array<double>^% c,
                           [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax);
  };
}