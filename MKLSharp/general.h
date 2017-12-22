#pragma once

using namespace System;
using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Lapack {
  public:
    static __int64 sgetrf(LapackLayout Layout, __int64 m, __int64 n,
                          array<float>^ a, __int64 lda, [Out]array<__int64>^% ipiv);
    static __int64 dgetrf(LapackLayout Layout, __int64 m, __int64 n,
                          array<double>^ a, __int64 lda, [Out]array<__int64>^% ipiv);

    static __int64 sgetrs(LapackLayout Layout, char trans,
                          __int64 n, __int64 nrhs, array<float>^ a, __int64 lda,
                          array<__int64>^ ipiv, array<float>^ b, __int64 ldb);
    static __int64 dgetrs(LapackLayout Layout, char trans,
                          __int64 n, __int64 nrhs, array<double>^ a, __int64 lda,
                          array<__int64>^ ipiv, array<double>^ b, __int64 ldb);
  };
}