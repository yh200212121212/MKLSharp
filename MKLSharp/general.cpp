#include "stdafx.h"

#include "general.h"

namespace MKLSharp {
  __int64 Lapack::sgetrf(LapackLayout Layout, __int64 m, __int64 n,
                         array<float>^ a, __int64 lda, [Out]array<__int64>^% ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    ipiv = gcnew array<__int64>((int)Math::Max(1LL, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_sgetrf((int)Layout, m, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dgetrf(LapackLayout Layout, __int64 m, __int64 n,
                         array<double>^ a, __int64 lda, [Out]array<__int64>^% ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    ipiv = gcnew array<__int64>((int)Math::Max(1LL, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dgetrf((int)Layout, m, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::sgetrs(LapackLayout Layout, char trans,
                         __int64 n, __int64 nrhs, array<float>^ a, __int64 lda,
                         array<__int64>^ ipiv, array<float>^ b, __int64 ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_sgetrs((int)Layout, trans, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dgetrs(LapackLayout Layout, char trans,
                         __int64 n, __int64 nrhs, array<double>^ a, __int64 lda,
                         array<__int64>^ ipiv, array<double>^ b, __int64 ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dgetrs((int)Layout, trans, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
}