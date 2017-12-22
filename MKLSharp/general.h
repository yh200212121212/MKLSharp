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
  };
}