#include "stdafx.h"

#include "general.h"

namespace MKLSharp {
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

  __int64 Lapack::sgetrs(LapackLayout Layout, char trans,
                         int n, int nrhs, array<float>^ a, int lda,
                         array<__int64>^ ipiv, array<float>^ b, int ldb) {
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
                         int n, int nrhs, array<double>^ a, int lda,
                         array<__int64>^ ipiv, array<double>^ b, int ldb) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<__int64> ptr_i = &ipiv[0];
    pin_ptr<double> ptr_b = &b[0];
    auto res = LAPACKE_dgetrs((int)Layout, trans, n, nrhs, ptr_a, lda, ptr_i, ptr_b, ldb);
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

  __int64 Lapack::sgecon(LapackLayout Layout, char norm,
                         int n, array<float>^ a, int lda,
                         float aNorm, [Out]float% rCond) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_rc = &rCond;
    auto res = LAPACKE_sgecon((int)Layout, norm, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }
  __int64 Lapack::dgecon(LapackLayout Layout, char norm,
                         int n, array<double>^ a, int lda,
                         double aNorm, [Out]double% rCond) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_rc = &rCond;
    auto res = LAPACKE_dgecon((int)Layout, norm, n, ptr_a, lda, aNorm, ptr_rc);
    ptr_a = nullptr;
    ptr_rc = nullptr;
    return res;
  }

  __int64 Lapack::sgerfs(LapackLayout Layout, char Trans,
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
    auto res = LAPACKE_sgerfs((int)Layout, Trans, n, nrhs, ptr_a, lda,
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
  __int64 Lapack::dgerfs(LapackLayout Layout, char Trans,
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
    auto res = LAPACKE_dgerfs((int)Layout, Trans, n, nrhs, ptr_a, lda,
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

  __int64 Lapack::sgerfsx(LapackLayout Layout, char trans, char equed,
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
    auto res = LAPACKE_sgerfsx((int)Layout, trans, equed,
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
  __int64 Lapack::dgerfsx(LapackLayout Layout, char trans, char equed,
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
    auto res = LAPACKE_dgerfsx((int)Layout, trans, equed,
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
}