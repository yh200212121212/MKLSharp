#include "stdafx.h"

#include "Lapack.h"

namespace MKLSharp {
  #pragma region general
  __int64 Lapack::sgetrf(LapackLayout Layout, int m, int n,
                         array<float>^ a, int lda, [Out]array<__int64>^% ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_sgetrf((int)Layout, m, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dgetrf(LapackLayout Layout, int m, int n,
                         array<double>^ a, int lda, [Out]array<__int64>^% ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dgetrf((int)Layout, m, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::sgetrs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<__int64>^ ipiv, array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_sgetrs((int)Layout, (char)Trans, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dgetrs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int nrhs, array<double>^ a, int lda,
                         array<__int64>^ ipiv, array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dgetrs((int)Layout, (char)Trans, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sgeequ(LapackLayout Layout, int m, int n,
                         array<float>^ a, int lda,
                         [Out]array<float>^% r, [Out]array<float>^% c,
                         [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &a[0];
    r = gcnew array<float>(m);
    pin_ptr<float> ptr_r = &r[0];
    c = gcnew array<float>(n);
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_rc = &rowCnd;
    pin_ptr<float> ptr_cc = &colCnd;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_sgeequ((int)Layout, m, n, ptr_a, lda,
                              ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dgeequ(LapackLayout Layout, int m, int n,
                         array<double>^ a, int lda,
                         [Out]array<double>^% r, [Out]array<double>^% c,
                         [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &a[0];
    r = gcnew array<double>(m);
    pin_ptr<double> ptr_r = &r[0];
    c = gcnew array<double>(n);
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_rc = &rowCnd;
    pin_ptr<double> ptr_cc = &colCnd;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dgeequ((int)Layout, m, n, ptr_a, lda,
                              ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::sgeequb(LapackLayout Layout, int m, int n,
                          array<float>^ a, int lda,
                          [Out]array<float>^% r, [Out]array<float>^% c,
                          [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &a[0];
    r = gcnew array<float>(m);
    pin_ptr<float> ptr_r = &r[0];
    c = gcnew array<float>(n);
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_rc = &rowCnd;
    pin_ptr<float> ptr_cc = &colCnd;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_sgeequb((int)Layout, m, n, ptr_a, lda,
                               ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dgeequb(LapackLayout Layout, int m, int n,
                          array<double>^ a, int lda,
                          [Out]array<double>^% r, [Out]array<double>^% c,
                          [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &a[0];
    r = gcnew array<double>(m);
    pin_ptr<double> ptr_r = &r[0];
    c = gcnew array<double>(n);
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_rc = &rowCnd;
    pin_ptr<double> ptr_cc = &colCnd;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dgeequb((int)Layout, m, n, ptr_a, lda,
                               ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::sgecon(LapackLayout Layout, LapackNorm Norm,
                         int n, array<float>^ a, int lda,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sgecon((int)Layout, (char)Norm, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dgecon(LapackLayout Layout, LapackNorm Norm,
                         int n, array<double>^ a, int lda,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dgecon((int)Layout, (char)Norm, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sgerfs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ af, int ldaf, array<__int64>^ ipiv,
                         array<float>^ b, int ldb, array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_sgerfs((int)Layout, (char)Trans, n, nrhs, ptr_a, lda,
                              ptr_af, ldaf, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dgerfs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ af, int ldaf, array<__int64>^ ipiv,
                         array<double>^ b, int ldb, array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dgerfs((int)Layout, (char)Trans, n, nrhs, ptr_a, lda,
                              ptr_af, ldaf, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::sgerfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf, array<__int64>^ ipiv,
                          array<float>^ r, array<float>^ c,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]float% rCond, [Out]array<float>^% bErr,
                          int nErrBnds, [Out]array<float>^% errBndsNorm,
                          [Out]array<float>^% errBndsComp,
                          int nParams, array<float>^ params) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_r = &r[0];
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_rc = &rCond;
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebc = &errBndsComp[0];
    pin_ptr<float> ptr_p = &params[0];
    auto res = LAPACKE_sgerfsx((int)Layout, (char)Trans, (char)Equed,
                               n, nrhs, ptr_a, lda,
                               ptr_af, ldaf, ptr_i,
                               ptr_r, ptr_c,
                               ptr_b, ldb,
                               ptr_x, ldx,
                               ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc,
                               nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }
  __int64 Lapack::dgerfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf, array<__int64>^ ipiv,
                          array<double>^ r, array<double>^ c,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]double% rCond, [Out]array<double>^% bErr,
                          int nErrBnds, [Out]array<double>^% errBndsNorm,
                          [Out]array<double>^% errBndsComp,
                          int nParams, array<double>^ params) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_r = &r[0];
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_rc = &rCond;
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebc = &errBndsComp[0];
    pin_ptr<double> ptr_p = &params[0];
    auto res = LAPACKE_dgerfsx((int)Layout, (char)Trans, (char)Equed,
                               n, nrhs, ptr_a, lda,
                               ptr_af, ldaf, ptr_i,
                               ptr_r, ptr_c,
                               ptr_b, ldb,
                               ptr_x, ldx,
                               ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc,
                               nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }

  __int64 Lapack::sgetri(LapackLayout Layout,
                         int n, array<float>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_sgetri((int)Layout, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dgetri(LapackLayout Layout,
                         int n, array<double>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dgetri((int)Layout, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region general band
  __int64 Lapack::sgbtrf(LapackLayout Layout, int m, int n, int kl, int ku,
                         array<float>^ ab, int ldab, [Out]array<__int64>^% ipiv) {
    pin_ptr<float> ptr_a = &ab[0];
    ipiv = gcnew array<__int64>(Math::Max(1, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_sgbtrf((int)Layout, m, n, kl, ku, ptr_a, ldab, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dgbtrf(LapackLayout Layout, int m, int n, int kl, int ku,
                         array<double>^ ab, int ldab, [Out]array<__int64>^% ipiv) {
    pin_ptr<double> ptr_a = &ab[0];
    ipiv = gcnew array<__int64>(Math::Max(1, Math::Min(m, n)));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dgbtrf((int)Layout, m, n, kl, ku, ptr_a, ldab, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::sgbtrs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int kl, int ku, int nrhs,
                         array<float>^ ab, int ldab, array<__int64>^ ipiv,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_sgbtrs((int)Layout, (char)Trans,
                              n, kl, ku, nrhs,
                              ptr_a, ldab, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dgbtrs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int kl, int ku, int nrhs,
                         array<double>^ ab, int ldab, array<__int64>^ ipiv,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dgbtrs((int)Layout, (char)Trans,
                              n, kl, ku, nrhs,
                              ptr_a, ldab, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sgbequ(LapackLayout Layout, int m, int n, int kl, int ku,
                         array<float>^ ab, int ldab,
                         [Out]array<float>^% r, [Out]array<float>^% c,
                         [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &ab[0];
    r = gcnew array<float>(m);
    pin_ptr<float> ptr_r = &r[0];
    c = gcnew array<float>(n);
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_rc = &rowCnd;
    pin_ptr<float> ptr_cc = &colCnd;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_sgbequ((int)Layout, m, n, kl, ku,
                              ptr_a, ldab, ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dgbequ(LapackLayout Layout, int m, int n, int kl, int ku,
                         array<double>^ ab, int ldab,
                         [Out]array<double>^% r, [Out]array<double>^% c,
                         [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &ab[0];
    r = gcnew array<double>(m);
    pin_ptr<double> ptr_r = &r[0];
    c = gcnew array<double>(n);
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_rc = &rowCnd;
    pin_ptr<double> ptr_cc = &colCnd;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dgbequ((int)Layout, m, n, kl, ku,
                              ptr_a, ldab, ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::sgbequb(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<float>^ ab, int ldab,
                          [Out]array<float>^% r, [Out]array<float>^% c,
                          [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &ab[0];
    r = gcnew array<float>(m);
    pin_ptr<float> ptr_r = &r[0];
    c = gcnew array<float>(n);
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_rc = &rowCnd;
    pin_ptr<float> ptr_cc = &colCnd;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_sgbequb((int)Layout, m, n, kl, ku,
                               ptr_a, ldab, ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dgbequb(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<double>^ ab, int ldab,
                          [Out]array<double>^% r, [Out]array<double>^% c,
                          [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &ab[0];
    r = gcnew array<double>(m);
    pin_ptr<double> ptr_r = &r[0];
    c = gcnew array<double>(n);
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_rc = &rowCnd;
    pin_ptr<double> ptr_cc = &colCnd;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dgbequb((int)Layout, m, n, kl, ku,
                               ptr_a, ldab, ptr_r, ptr_c, ptr_rc, ptr_cc, ptr_am);
    ptr_a = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_rc = nullptr;
    ptr_cc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::sgbcon(LapackLayout Layout, LapackNorm Norm, int n, int kl, int ku,
                         array<float>^ ab, int ldab, array<__int64>^ ipiv,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sgbcon((int)Layout, (char)Norm, n, kl, ku,
                              ptr_a, ldab, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dgbcon(LapackLayout Layout, LapackNorm Norm, int n, int kl, int ku,
                         array<double>^ ab, int ldab, array<__int64>^ ipiv,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dgbcon((int)Layout, (char)Norm, n, kl, ku,
                              ptr_a, ldab, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sgbrfs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int kl, int ku, int nrhs, array<float>^ ab, int ldab,
                         array<float>^ afb, int ldafb, array<__int64>^ ipiv,
                         array<float>^ b, int ldb, 
                         array<float>^ x, int ldx, 
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<float> ptr_af = &afb[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_sgbrfs((int)Layout, (char)Trans, n, kl, ku, nrhs, ptr_a, ldab,
                              ptr_af, ldafb, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_i = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dgbrfs(LapackLayout Layout, LapackTranspose Trans,
                         int n, int kl, int ku, int nrhs, array<double>^ ab, int ldab,
                         array<double>^ afb, int ldafb, array<__int64>^ ipiv,
                         array<double>^ b, int ldb, 
                         array<double>^ x, int ldx, 
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<double> ptr_af = &afb[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dgbrfs((int)Layout, (char)Trans, n, kl, ku, nrhs, ptr_a, ldab,
                              ptr_af, ldafb, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_i = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::sgbrfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                          int n, int kl, int ku, int nrhs, array<float>^ ab, int ldab,
                          array<float>^ afb, int ldafb, array<__int64>^ ipiv,
                          array<float>^ r, array<float>^ c,
                          array<float>^ b, int ldb, 
                          array<float>^ x, int ldx, 
                          [Out]float% rCond, [Out]array<float>^% bErr,
                          int nErrBnds, [Out]array<float>^% errBndsNorm,
                          [Out]array<float>^% errBndsComp, int nParams, array<float>^ params) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<float> ptr_af = &afb[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_r = &r[0];
    pin_ptr<float> ptr_c = &c[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_rc = &rCond;
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebc = &errBndsComp[0];
    pin_ptr<float> ptr_p = &params[0];
    auto res = LAPACKE_sgbrfsx((int)Layout, (char)Trans, (char)Equed, 
                               n, kl, ku, nrhs, ptr_a, ldab,
                               ptr_af, ldafb, ptr_i,
                               ptr_r, ptr_c,
                               ptr_b, ldb,
                               ptr_x, ldx,
                               ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc,
                               nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }
  __int64 Lapack::dgbrfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                          int n, int kl, int ku, int nrhs, array<double>^ ab, int ldab,
                          array<double>^ afb, int ldafb, array<__int64>^ ipiv,
                          array<double>^ r, array<double>^ c,
                          array<double>^ b, int ldb, 
                          array<double>^ x, int ldx, 
                          [Out]double% rCond, [Out]array<double>^% bErr,
                          int nErrBnds, [Out]array<double>^% errBndsNorm,
                          [Out]array<double>^% errBndsComp, int nParams, array<double>^ params) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<double> ptr_af = &afb[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_r = &r[0];
    pin_ptr<double> ptr_c = &c[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_rc = &rCond;
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebc = &errBndsComp[0];
    pin_ptr<double> ptr_p = &params[0];
    auto res = LAPACKE_dgbrfsx((int)Layout, (char)Trans, (char)Equed, 
                               n, kl, ku, nrhs, ptr_a, ldab,
                               ptr_af, ldafb, ptr_i,
                               ptr_r, ptr_c,
                               ptr_b, ldb,
                               ptr_x, ldx,
                               ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc,
                               nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_r = nullptr;
    ptr_c = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region general tridiagonal
  __int64 Lapack::sgttrf(int n, array<float>^ dl, array<float>^ d, array<float>^ du,
                         [Out]array<float>^% du2, [Out]array<__int64>^% ipiv) {
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    du2 = gcnew array<float>(n - 2);
    pin_ptr<float> ptr_u2 = &du2[0];
    ipiv = gcnew array<__int64>(n);
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_sgttrf(n, ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dgttrf(int n, array<double>^ dl, array<double>^ d, array<double>^ du,
                         [Out]array<double>^% du2, [Out]array<__int64>^% ipiv) {
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    du2 = gcnew array<double>(n - 2);
    pin_ptr<double> ptr_u2 = &du2[0];
    ipiv = gcnew array<__int64>(n);
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dgttrf(n, ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::sgttrs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                         array<float>^ dl, array<float>^ d, array<float>^ du,
                         array<float>^ du2, array<__int64>^ ipiv,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    pin_ptr<float> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_sgttrs((int)Layout, (char)Trans, n, nrhs,
                              ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i, ptr_b, ldb);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dgttrs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                         array<double>^ dl, array<double>^ d, array<double>^ du,
                         array<double>^ du2, array<__int64>^ ipiv,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    pin_ptr<double> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dgttrs((int)Layout, (char)Trans, n, nrhs,
                              ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i, ptr_b, ldb);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sgtcon(LapackNorm Norm, int n,
                         array<float>^ dl, array<float>^ d, array<float>^ du,
                         array<float>^ du2, array<__int64>^ ipiv,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    pin_ptr<float> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sgtcon((char)Norm, n, ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i, aNorm, ptr_rc);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dgtcon(LapackNorm Norm, int n,
                         array<double>^ dl, array<double>^ d, array<double>^ du,
                         array<double>^ du2, array<__int64>^ ipiv,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    pin_ptr<double> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dgtcon((char)Norm, n, ptr_l, ptr_d, ptr_u, ptr_u2, ptr_i, aNorm, ptr_rc);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sgtrfs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                         array<float>^ dl, array<float>^ d, array<float>^ du,
                         array<float>^ dlf, array<float>^ df, array<float>^ duf,
                         array<float>^ du2, array<__int64>^ ipiv,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    pin_ptr<float> ptr_lf = &dlf[0];
    pin_ptr<float> ptr_df = &df[0];
    pin_ptr<float> ptr_uf = &duf[0];
    pin_ptr<float> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_sgtrfs((int)Layout, (char)Trans, n, nrhs,
                              ptr_l, ptr_d, ptr_u, ptr_lf, ptr_df, ptr_uf,
                              ptr_u2, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_lf = nullptr;
    ptr_df = nullptr;
    ptr_uf = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dgtrfs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                         array<double>^ dl, array<double>^ d, array<double>^ du,
                         array<double>^ dlf, array<double>^ df, array<double>^ duf,
                         array<double>^ du2, array<__int64>^ ipiv,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    pin_ptr<double> ptr_lf = &dlf[0];
    pin_ptr<double> ptr_df = &df[0];
    pin_ptr<double> ptr_uf = &duf[0];
    pin_ptr<double> ptr_u2 = &du2[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dgtrfs((int)Layout, (char)Trans, n, nrhs,
                              ptr_l, ptr_d, ptr_u, ptr_lf, ptr_df, ptr_uf,
                              ptr_u2, ptr_i, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_lf = nullptr;
    ptr_df = nullptr;
    ptr_uf = nullptr;
    ptr_u2 = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region diagonally dominant tridiagonal
  __int64 Lapack::sdttrfb(int n, array<float>^ dl, array<float>^ d, array<float>^ du) {
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    __int64 info, ln = n;
    SDTTRFB(&ln, ptr_l, ptr_d, ptr_u, &info);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    return info;
  }
  __int64 Lapack::ddttrfb(int n, array<double>^ dl, array<double>^ d, array<double>^ du) {
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    __int64 info, ln = n;
    DDTTRFB(&ln, ptr_l, ptr_d, ptr_u, &info);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    return info;
  }

  __int64 Lapack::sdttrsb(LapackTranspose Trans, int n, int nrhs,
                          array<float>^ dl, array<float>^ d, array<float>^ du,
                          array<float>^ b, int ldb) {
    char trans = (char)Trans;
    __int64 ln = n;
    __int64 lnrhs = nrhs;
    pin_ptr<float> ptr_l = &dl[0];
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_u = &du[0];
    pin_ptr<float> ptr_b = &b[0];
    __int64 lldb = ldb;
    __int64 info;
    SDTTRSB(&trans, &ln, &lnrhs, ptr_l, ptr_d, ptr_u, ptr_b, &lldb, &info);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_b = nullptr;
    return info;
  }
  __int64 Lapack::ddttrsb(LapackTranspose Trans, int n, int nrhs,
                          array<double>^ dl, array<double>^ d, array<double>^ du,
                          array<double>^ b, int ldb) {
    char trans = (char)Trans;
    __int64 ln = n;
    __int64 lnrhs = nrhs;
    pin_ptr<double> ptr_l = &dl[0];
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_u = &du[0];
    pin_ptr<double> ptr_b = &b[0];
    __int64 lldb = ldb;
    __int64 info;
    DDTTRSB(&trans, &ln, &lnrhs, ptr_l, ptr_d, ptr_u, ptr_b, &lldb, &info);
    ptr_l = nullptr;
    ptr_d = nullptr;
    ptr_u = nullptr;
    ptr_b = nullptr;
    return info;
  }
  #pragma endregion
  #pragma region symmetric positive-definite
  __int64 Lapack::spotrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda) {
    pin_ptr<float> ptr_a = &a[0];
    auto res = LAPACKE_spotrf((int)Layout, (char)UpLo, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpotrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda) {
    pin_ptr<double> ptr_a = &a[0];
    auto res = LAPACKE_dpotrf((int)Layout, (char)UpLo, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }

  __int64 Lapack::spotrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_spotrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dpotrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dpotrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::spoequ(LapackLayout Layout, int n,
                         array<float>^ a, int lda,
                         [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &a[0];
    s = gcnew array<float>(n);
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_sc = &sCond;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_spoequ((int)Layout, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dpoequ(LapackLayout Layout, int n,
                         array<double>^ a, int lda,
                         [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &a[0];
    s = gcnew array<double>(n);
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_sc = &sCond;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dpoequ((int)Layout, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::spoequb(LapackLayout Layout, int n,
                          array<float>^ a, int lda,
                          [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &a[0];
    s = gcnew array<float>(n);
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_sc = &sCond;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_spoequb((int)Layout, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dpoequb(LapackLayout Layout, int n,
                          array<double>^ a, int lda,
                          [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &a[0];
    s = gcnew array<double>(n);
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_sc = &sCond;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dpoequb((int)Layout, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::spocon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_spocon((int)Layout, (char)UpLo, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dpocon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dpocon((int)Layout, (char)UpLo, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sporfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ af, int ldaf,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_sporfs((int)Layout, (char)UpLo,
                              n, nrhs, ptr_a, lda, ptr_af, ldaf,
                              ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dporfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ af, int ldaf,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dporfs((int)Layout, (char)UpLo,
                              n, nrhs, ptr_a, lda, ptr_af, ldaf,
                              ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::sporfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf, array<float>^ s,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]float% rCond, [Out]array<float>^% bErr,
                          int nErrBnds, [Out]array<float>^% errBndsNorm,
                          [Out]array<float>^% errBndsComp,
                          int nParams, array<float>^ params) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_rc = &rCond;
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebc = &errBndsComp[0];
    pin_ptr<float> ptr_p = &params[0];
    auto res = LAPACKE_sporfsx((int)Layout, (char)UpLo, (char)Equed,
                               n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_s,
                               ptr_b, ldb, ptr_x, ldx, ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc, nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_s = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }
  __int64 Lapack::dporfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf, array<double>^ s,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]double% rCond, [Out]array<double>^% bErr,
                          int nErrBnds, [Out]array<double>^% errBndsNorm,
                          [Out]array<double>^% errBndsComp,
                          int nParams, array<double>^ params) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_rc = &rCond;
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebc = &errBndsComp[0];
    pin_ptr<double> ptr_p = &params[0];
    auto res = LAPACKE_dporfsx((int)Layout, (char)UpLo, (char)Equed,
                               n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_s,
                               ptr_b, ldb, ptr_x, ldx, ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc, nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_s = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }

  __int64 Lapack::spotri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda) {
    pin_ptr<float> ptr_a = &a[0];
    auto res = LAPACKE_spotri((int)Layout, (char)UpLo, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpotri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda) {
    pin_ptr<double> ptr_a = &a[0];
    auto res = LAPACKE_dpotri((int)Layout, (char)UpLo, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric positive-definite, packed storage
  __int64 Lapack::spptrf(LapackLayout Layout, LapackUpLo UpLo, int n, array<float>^ ap) {
    pin_ptr<float> ptr_a = &ap[0];
    auto res = LAPACKE_spptrf((int)Layout, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpptrf(LapackLayout Layout, LapackUpLo UpLo, int n, array<double>^ ap) {
    pin_ptr<double> ptr_a = &ap[0];
    auto res = LAPACKE_dpptrf((int)Layout, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }

  __int64 Lapack::spptrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ ap, array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_spptrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dpptrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ ap, array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dpptrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sppequ(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ ap,
                         [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &ap[0];
    s = gcnew array<float>(n);
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_sc = &sCond;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_sppequ((int)Layout, (char)UpLo, n, ptr_a, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dppequ(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ ap,
                         [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &ap[0];
    s = gcnew array<double>(n);
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_sc = &sCond;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dppequ((int)Layout, (char)UpLo, n, ptr_a, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::sppcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ ap, float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sppcon((int)Layout, (char)UpLo, n, ptr_a, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dppcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ ap, double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dppcon((int)Layout, (char)UpLo, n, ptr_a, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::spprfs(LapackLayout Layout, LapackUpLo UpLo, 
                         int n, int nrhs, array<float>^ ap, array<float>^ afp,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_af = &afp[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_spprfs((int)Layout, (char)UpLo,
                              n, nrhs, ptr_a, ptr_af, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dpprfs(LapackLayout Layout, LapackUpLo UpLo, 
                         int n, int nrhs, array<double>^ ap, array<double>^ afp,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_af = &afp[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dpprfs((int)Layout, (char)UpLo,
                              n, nrhs, ptr_a, ptr_af, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::spptri(LapackLayout Layout, LapackUpLo UpLo, int n, array<float>^ ap) {
    pin_ptr<float> ptr_a = &ap[0];
    auto res = LAPACKE_spptri((int)Layout, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpptri(LapackLayout Layout, LapackUpLo UpLo, int n, array<double>^ ap) {
    pin_ptr<double> ptr_a = &ap[0];
    auto res = LAPACKE_dpptri((int)Layout, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric positive-definite, RFP storage
  __int64 Lapack::spftrf(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, array<float>^ a) {
    pin_ptr<float> ptr_a = &a[0];
    auto res = LAPACKE_spftrf((int)Layout, (char)Trans, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpftrf(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, array<double>^ a) {
    pin_ptr<double> ptr_a = &a[0];
    auto res = LAPACKE_dpftrf((int)Layout, (char)Trans, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }

  __int64 Lapack::spftrs(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ a, array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_spftrs((int)Layout, (char)Trans, (char)UpLo, n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dpftrs(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ a, array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dpftrs((int)Layout, (char)Trans, (char)UpLo, n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::spftri(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, array<float>^ a) {
    pin_ptr<float> ptr_a = &a[0];
    auto res = LAPACKE_spftri((int)Layout, (char)Trans, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpftri(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                         int n, array<double>^ a) {
    pin_ptr<double> ptr_a = &a[0];
    auto res = LAPACKE_dpftri((int)Layout, (char)Trans, (char)UpLo, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric positive-definite, band
  __int64 Lapack::spbtrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<float>^ ab, int ldab) {
    pin_ptr<float> ptr_a = &ab[0];
    auto res = LAPACKE_spbtrf((int)Layout, (char)UpLo, n, kd, ptr_a, ldab);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dpbtrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<double>^ ab, int ldab) {
    pin_ptr<double> ptr_a = &ab[0];
    auto res = LAPACKE_dpbtrf((int)Layout, (char)UpLo, n, kd, ptr_a, ldab);
    ptr_a = nullptr;
    return res;
  }

  __int64 Lapack::spbtrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, int nrhs, array<float>^ ab, int ldab,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_spbtrs((int)Layout, (char)UpLo, n, kd, nrhs, ptr_a, ldab, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dpbtrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, int nrhs, array<double>^ ab, int ldab,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dpbtrs((int)Layout, (char)UpLo, n, kd, nrhs, ptr_a, ldab, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::spbequ(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<float>^ ab, int ldab,
                         [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &ab[0];
    s = gcnew array<float>(n);
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_sc = &sCond;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_spbequ((int)Layout, (char)UpLo, n, kd, ptr_a, ldab,
                              ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dpbequ(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<double>^ ab, int ldab,
                         [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &ab[0];
    s = gcnew array<double>(n);
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_sc = &sCond;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dpbequ((int)Layout, (char)UpLo, n, kd, ptr_a, ldab,
                              ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::spbcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<float>^ ab, int ldab,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_spbcon((int)Layout, (char)UpLo, n, kd, ptr_a, ldab,
                              aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dpbcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, array<double>^ ab, int ldab,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dpbcon((int)Layout, (char)UpLo, n, kd, ptr_a, ldab,
                              aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::spbrfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, int nrhs, array<float>^ ab, int ldab,
                         array<float>^ afb, int ldafb, 
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &ab[0];
    pin_ptr<float> ptr_af = &afb[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_spbrfs((int)Layout, (char)UpLo, n, kd, nrhs, ptr_a, ldab, ptr_af, ldafb,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dpbrfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int kd, int nrhs, array<double>^ ab, int ldab,
                         array<double>^ afb, int ldafb, 
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &ab[0];
    pin_ptr<double> ptr_af = &afb[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dpbrfs((int)Layout, (char)UpLo, n, kd, nrhs, ptr_a, ldab, ptr_af, ldafb,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric positive-definite, tridiagonal
  __int64 Lapack::spttrf(int n, array<float>^ d, array<float>^ e) {
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_e = &e[0];
    auto res = LAPACKE_spttrf(n, ptr_d, ptr_e);
    ptr_d = nullptr;
    ptr_e = nullptr;
    return res;
  }
  __int64 Lapack::dpttrf(int n, array<double>^ d, array<double>^ e) {
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_e = &e[0];
    auto res = LAPACKE_dpttrf(n, ptr_d, ptr_e);
    ptr_d = nullptr;
    ptr_e = nullptr;
    return res;
  }

  __int64 Lapack::spttrs(LapackLayout Layout, int n, int nrhs,
                         array<float>^ d, array<float>^ e,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_e = &e[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_spttrs((int)Layout, n, nrhs, ptr_d, ptr_e, ptr_b, ldb);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dpttrs(LapackLayout Layout, int n, int nrhs,
                         array<double>^ d, array<double>^ e,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_e = &e[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dpttrs((int)Layout, n, nrhs, ptr_d, ptr_e, ptr_b, ldb);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sptcon(int n, array<float>^ d, array<float>^ e,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_e = &e[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sptcon(n, ptr_d, ptr_e, aNorm, ptr_rc);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dptcon(int n, array<double>^ d, array<double>^ e,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_e = &e[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dptcon(n, ptr_d, ptr_e, aNorm, ptr_rc);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sptrfs(LapackLayout Layout, int n, int nrhs,
                         array<float>^ d, array<float>^ e,
                         array<float>^ df, array<float>^ ef,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_d = &d[0];
    pin_ptr<float> ptr_e = &e[0];
    pin_ptr<float> ptr_df = &df[0];
    pin_ptr<float> ptr_ef = &ef[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_sptrfs((int)Layout, n, nrhs, ptr_d, ptr_e, ptr_df, ptr_ef,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_df = nullptr;
    ptr_ef = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dptrfs(LapackLayout Layout, int n, int nrhs,
                         array<double>^ d, array<double>^ e,
                         array<double>^ df, array<double>^ ef,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_d = &d[0];
    pin_ptr<double> ptr_e = &e[0];
    pin_ptr<double> ptr_df = &df[0];
    pin_ptr<double> ptr_ef = &ef[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dptrfs((int)Layout, n, nrhs, ptr_d, ptr_e, ptr_df, ptr_ef,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_d = nullptr;
    ptr_e = nullptr;
    ptr_df = nullptr;
    ptr_ef = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric indefinite
  __int64 Lapack::ssytrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda, [Out]array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, n));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssytrf((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsytrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda, [Out]array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, n));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsytrf((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::ssytrf_aa(LapackLayout Layout, LapackUpLo UpLo,
                            int n, array<float>^ a, int lda, [Out]array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, n));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssytrf_aa((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsytrf_aa(LapackLayout Layout, LapackUpLo UpLo,
                            int n, array<double>^ a, int lda, [Out]array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    ipiv = gcnew array<__int64>(Math::Max(1, n));
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsytrf_aa((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::ssytrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_ssytrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dsytrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dsytrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::ssytrs2(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                          array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_ssytrs2((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dsytrs2(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                          array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dsytrs2((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
 
  __int64 Lapack::ssytrs_aa(LapackLayout Layout, LapackUpLo UpLo,
                            int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                            array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_ssytrs_aa((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dsytrs_aa(LapackLayout Layout, LapackUpLo UpLo,
                            int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                            array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dsytrs_aa((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::ssyequb(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda,
                          [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax) {
    pin_ptr<float> ptr_a = &a[0];
    s = gcnew array<float>(n);
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_sc = &sCond;
    pin_ptr<float> ptr_am = &aMax;
    auto res = LAPACKE_ssyequb((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }
  __int64 Lapack::dsyequb(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda,
                          [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax) {
    pin_ptr<double> ptr_a = &a[0];
    s = gcnew array<double>(n);
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_sc = &sCond;
    pin_ptr<double> ptr_am = &aMax;
    auto res = LAPACKE_dsyequb((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_s, ptr_sc, ptr_am);
    ptr_a = nullptr;
    ptr_s = nullptr;
    ptr_sc = nullptr;
    ptr_am = nullptr;
    return res;
  }

  __int64 Lapack::ssycon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda, array<__int64>^ ipiv,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_ssycon((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dsycon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda, array<__int64>^ ipiv,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dsycon((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::ssyrfs(LapackLayout Layout, LapackUpLo UpLo, 
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ af, int ldaf, array<__int64>^ ipiv,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx, 
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_ssyrfs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_i,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dsyrfs(LapackLayout Layout, LapackUpLo UpLo, 
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ af, int ldaf, array<__int64>^ ipiv,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx, 
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dsyrfs((int)Layout, (char)UpLo, n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_i,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::ssyrfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf, array<__int64>^ ipiv,
                          array<float>^ s, array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]float% rCond, [Out]array<float>^% bErr, 
                          int nErrBnds, [Out]array<float>^% errBndsNorm, 
                          [Out]array<float>^% errBndsComp,
                          int nParams, array<float>^ params) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_s = &s[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_rc = &rCond;
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<float>(nrhs * nErrBnds);
    pin_ptr<float> ptr_ebc = &errBndsComp[0];
    pin_ptr<float> ptr_p = &params[0];
    auto res = LAPACKE_ssyrfsx((int)Layout, (char)UpLo, (char)Equed,
                               n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_i,
                               ptr_s, ptr_b, ldb, ptr_x, ldx, ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc, nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_s = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }
  __int64 Lapack::dsyrfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf, array<__int64>^ ipiv,
                          array<double>^ s, array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]double% rCond, [Out]array<double>^% bErr, 
                          int nErrBnds, [Out]array<double>^% errBndsNorm, 
                          [Out]array<double>^% errBndsComp,
                          int nParams, array<double>^ params) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_af = &af[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_s = &s[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_rc = &rCond;
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    errBndsNorm = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebn = &errBndsNorm[0];
    errBndsComp = gcnew array<double>(nrhs * nErrBnds);
    pin_ptr<double> ptr_ebc = &errBndsComp[0];
    pin_ptr<double> ptr_p = &params[0];
    auto res = LAPACKE_dsyrfsx((int)Layout, (char)UpLo, (char)Equed,
                               n, nrhs, ptr_a, lda, ptr_af, ldaf, ptr_i,
                               ptr_s, ptr_b, ldb, ptr_x, ldx, ptr_rc, ptr_be,
                               nErrBnds, ptr_ebn, ptr_ebc, nParams, ptr_p);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_s = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_rc = nullptr;
    ptr_be = nullptr;
    ptr_ebn = nullptr;
    ptr_ebc = nullptr;
    ptr_p = nullptr;
    return res;
  }

  __int64 Lapack::ssytri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssytri((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsytri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsytri((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::ssytri2(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssytri2((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsytri2(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda, array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsytri2((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  __int64 Lapack::ssytri2x(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<float>^ a, int lda, array<__int64>^ ipiv, int nb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssytri2x((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i, nb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsytri2x(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<double>^ a, int lda, array<__int64>^ ipiv, int nb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsytri2x((int)Layout, (char)UpLo, n, ptr_a, lda, ptr_i, nb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region symmetric indefinite, packed storage
  __int64 Lapack::ssptrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ ap, [Out]array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &ap[0];
    ipiv = gcnew array<__int64>(n);
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssptrf((int)Layout, (char)UpLo, n, ptr_a, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsptrf(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ ap, [Out]array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &ap[0];
    ipiv = gcnew array<__int64>(n);
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsptrf((int)Layout, (char)UpLo, n, ptr_a, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }

  void Lapack::sspffrt2(int n, int nColumn, array<float>^ ap) {
    pin_ptr<float> ptr_a = &ap[0];
    __int64 ln = n;
    __int64 lnc = nColumn;
    auto work = new float[n];
    auto work2 = new float[n];
    mkl_sspffrt2(ptr_a, &ln, &lnc, work, work2);
    delete work;
    delete work2;
    ptr_a = nullptr;
  }
  void Lapack::dspffrt2(int n, int nColumn, array<double>^ ap) {
    pin_ptr<double> ptr_a = &ap[0];
    __int64 ln = n;
    __int64 lnc = nColumn;
    auto work = new double[n];
    auto work2 = new double[n];
    mkl_dspffrt2(ptr_a, &ln, &lnc, work, work2);
    delete work;
    delete work2;
    ptr_a = nullptr;
  }

  void Lapack::sspffrtx(int n, int nColumn, array<float>^ ap) {
    pin_ptr<float> ptr_a = &ap[0];
    __int64 ln = n;
    __int64 lnc = nColumn;
    auto work = new float[n];
    auto work2 = new float[n];
    mkl_sspffrtx(ptr_a, &ln, &lnc, work, work2);
    delete work;
    delete work2;
    ptr_a = nullptr;
  }
  void Lapack::dspffrtx(int n, int nColumn, array<double>^ ap) {
    pin_ptr<double> ptr_a = &ap[0];
    __int64 ln = n;
    __int64 lnc = nColumn;
    auto work = new double[n];
    auto work2 = new double[n];
    mkl_dspffrtx(ptr_a, &ln, &lnc, work, work2);
    delete work;
    delete work2;
    ptr_a = nullptr;
  }

  __int64 Lapack::ssptrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ ap, array<__int64>^ ipiv,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_ssptrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dsptrs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ ap, array<__int64>^ ipiv,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dsptrs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_i, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::sspcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ ap, array<__int64>^ ipiv,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sspcon((int)Layout, (char)UpLo, n, ptr_a, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dspcon(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ ap, array<__int64>^ ipiv,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dspcon((int)Layout, (char)UpLo, n, ptr_a, ptr_i, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_i = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::ssprfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<float>^ ap,
                         array<float>^ afp, array<__int64>^ ipiv,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx,
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_af = &afp[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(Math::Max(1, nrhs));
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_ssprfs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_af, ptr_i,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dsprfs(LapackLayout Layout, LapackUpLo UpLo,
                         int n, int nrhs, array<double>^ ap,
                         array<double>^ afp, array<__int64>^ ipiv,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx,
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_af = &afp[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(Math::Max(1, nrhs));
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dsprfs((int)Layout, (char)UpLo, n, nrhs, ptr_a, ptr_af, ptr_i,
                              ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_af = nullptr;
    ptr_i = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::ssptri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<float>^ ap, array<__int64>^ ipiv) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_ssptri((int)Layout, (char)UpLo, n, ptr_a, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  __int64 Lapack::dsptri(LapackLayout Layout, LapackUpLo UpLo,
                         int n, array<double>^ ap, array<__int64>^ ipiv) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    auto res = LAPACKE_dsptri((int)Layout, (char)UpLo, n, ptr_a, ptr_i);
    ptr_a = nullptr;
    ptr_i = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region triangular
  __int64 Lapack::strtrs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_strtrs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, lda, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dtrtrs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag,
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dtrtrs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, lda, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::strcon(LapackLayout Layout, LapackNorm Norm,
                         LapackUpLo UpLo, LapackDiag Diag, 
                         int n, array<float>^ a, int lda, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_strcon((int)Layout, (char)Norm, (char)UpLo, (char)Diag,
                              n, ptr_a, lda, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dtrcon(LapackLayout Layout, LapackNorm Norm,
                         LapackUpLo UpLo, LapackDiag Diag, 
                         int n, array<double>^ a, int lda, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dtrcon((int)Layout, (char)Norm, (char)UpLo, (char)Diag,
                              n, ptr_a, lda, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::strrfs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag, 
                         int n, int nrhs, array<float>^ a, int lda,
                         array<float>^ b, int ldb,
                         array<float>^ x, int ldx, 
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(nrhs > 1 ? nrhs : 1);
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(nrhs > 1 ? nrhs : 1);
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_strrfs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, lda, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dtrrfs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag, 
                         int n, int nrhs, array<double>^ a, int lda,
                         array<double>^ b, int ldb,
                         array<double>^ x, int ldx, 
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(nrhs > 1 ? nrhs : 1);
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(nrhs > 1 ? nrhs : 1);
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dtrrfs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, lda, ptr_b, ldb, ptr_x, ldx,
                              ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::strtri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<float>^ a, int lda) {
    pin_ptr<float> ptr_a = &a[0];
    auto res = LAPACKE_strtri((int)Layout, (char)UpLo, (char)Diag, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dtrtri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<double>^ a, int lda) {
    pin_ptr<double> ptr_a = &a[0];
    auto res = LAPACKE_dtrtri((int)Layout, (char)UpLo, (char)Diag, n, ptr_a, lda);
    ptr_a = nullptr;
    return res;
  }
  #pragma endregion
  #pragma region triangular, packed storage
  __int64 Lapack::stptrs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag,
                         int n, int nrhs, array<float>^ ap,
                         array<float>^ b, int ldb) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_b = &b[0];
    auto res = LAPACKE_stptrs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }
  __int64 Lapack::dtptrs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag,
                         int n, int nrhs, array<double>^ ap,
                         array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dtptrs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, ptr_b, ldb);
    ptr_a = nullptr;
    ptr_b = nullptr;
    return res;
  }

  __int64 Lapack::stpcon(LapackLayout Layout, LapackNorm Norm,
                         LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<float>^ ap, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_stpcon((int)Layout, (char)Norm, (char)UpLo, (char)Diag,
                              n, ptr_a, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dtpcon(LapackLayout Layout, LapackNorm Norm,
                         LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<double>^ ap, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dtpcon((int)Layout, (char)Norm, (char)UpLo, (char)Diag,
                              n, ptr_a, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::stprfs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag, 
                         int n, int nrhs, array<float>^ ap,
                         array<float>^ b, int ldb, 
                         array<float>^ x, int ldx, 
                         [Out]array<float>^% fErr, [Out]array<float>^% bErr) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_x = &x[0];
    fErr = gcnew array<float>(nrhs > 1 ? nrhs : 1);
    pin_ptr<float> ptr_fe = &fErr[0];
    bErr = gcnew array<float>(nrhs > 1 ? nrhs : 1);
    pin_ptr<float> ptr_be = &bErr[0];
    auto res = LAPACKE_stprfs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }
  __int64 Lapack::dtprfs(LapackLayout Layout, LapackUpLo UpLo,
                         LapackTranspose Trans, LapackDiag Diag, 
                         int n, int nrhs, array<double>^ ap,
                         array<double>^ b, int ldb, 
                         array<double>^ x, int ldx, 
                         [Out]array<double>^% fErr, [Out]array<double>^% bErr) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_x = &x[0];
    fErr = gcnew array<double>(nrhs > 1 ? nrhs : 1);
    pin_ptr<double> ptr_fe = &fErr[0];
    bErr = gcnew array<double>(nrhs > 1 ? nrhs : 1);
    pin_ptr<double> ptr_be = &bErr[0];
    auto res = LAPACKE_dtprfs((int)Layout, (char)UpLo, (char)Trans, (char)Diag,
                              n, nrhs, ptr_a, ptr_b, ldb, ptr_x, ldx, ptr_fe, ptr_be);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_x = nullptr;
    ptr_fe = nullptr;
    ptr_be = nullptr;
    return res;
  }

  __int64 Lapack::stptri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<float>^ ap) {
    pin_ptr<float> ptr_a = &ap[0];
    auto res = LAPACKE_stptri((int)Layout, (char)UpLo, (char)Diag, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  __int64 Lapack::dtptri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                         int n, array<double>^ ap) {
    pin_ptr<double> ptr_a = &ap[0];
    auto res = LAPACKE_dtptri((int)Layout, (char)UpLo, (char)Diag, n, ptr_a);
    ptr_a = nullptr;
    return res;
  }
  #pragma endregion
}