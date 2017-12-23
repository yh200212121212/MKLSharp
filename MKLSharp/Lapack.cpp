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
  #pragma endregion
}