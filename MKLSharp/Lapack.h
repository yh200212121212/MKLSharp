#pragma once

using namespace System;
using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Lapack {
  public:
  #pragma region general
    static __int64 sgetrf(LapackLayout Layout, int m, int n,
                          array<float>^ a, int lda, [Out]array<__int64>^% ipiv);
    static __int64 dgetrf(LapackLayout Layout, int m, int n,
                          array<double>^ a, int lda, [Out]array<__int64>^% ipiv);

    static __int64 sgetrs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<__int64>^ ipiv, array<float>^ b, int ldb);
    static __int64 dgetrs(LapackLayout Layout, LapackTranspose Trans,
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

    static __int64 sgecon(LapackLayout Layout, LapackNorm Norm,
                          int n, array<float>^ a, int lda,
                          float aNorm, [Out]float% rCond);
    static __int64 dgecon(LapackLayout Layout, LapackNorm Norm,
                          int n, array<double>^ a, int lda,
                          double aNorm, [Out]double% rCond);

    static __int64 sgerfs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf, array<__int64>^ ipiv,
                          array<float>^ b, int ldb, array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dgerfs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf, array<__int64>^ ipiv,
                          array<double>^ b, int ldb, array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 sgerfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                           int n, int nrhs, array<float>^ a, int lda,
                           array<float>^ af, int ldaf, array<__int64>^ ipiv,
                           array<float>^ r, array<float>^ c,
                           array<float>^ b, int ldb,
                           array<float>^ x, int ldx,
                           [Out]float% rCond, [Out]array<float>^% bErr,
                           int nErrBnds, [Out]array<float>^% errBndsNorm,
                           [Out]array<float>^% errBndsComp,
                           int nParams, array<float>^ params);
    static __int64 dgerfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                           int n, int nrhs, array<double>^ a, int lda,
                           array<double>^ af, int ldaf, array<__int64>^ ipiv,
                           array<double>^ r, array<double>^ c,
                           array<double>^ b, int ldb,
                           array<double>^ x, int ldx,
                           [Out]double% rCond, [Out]array<double>^% bErr,
                           int nErrBnds, [Out]array<double>^% errBndsNorm,
                           [Out]array<double>^% errBndsComp,
                           int nParams, array<double>^ params);

    static __int64 sgetri(LapackLayout Layout,
                          int n, array<float>^ a, int lda, array<__int64>^ ipiv);
    static __int64 dgetri(LapackLayout Layout,
                          int n, array<double>^ a, int lda, array<__int64>^ ipiv);
  #pragma endregion
  #pragma region general band
    static __int64 sgbtrf(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<float>^ ab, int ldab, [Out]array<__int64>^% ipiv);
    static __int64 dgbtrf(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<double>^ ab, int ldab, [Out]array<__int64>^% ipiv);

    static __int64 sgbtrs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int kl, int ku, int nrhs,
                          array<float>^ ab, int ldab, array<__int64>^ ipiv,
                          array<float>^ b, int ldb);
    static __int64 dgbtrs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int kl, int ku, int nrhs,
                          array<double>^ ab, int ldab, array<__int64>^ ipiv,
                          array<double>^ b, int ldb);

    static __int64 sgbequ(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<float>^ ab, int ldab,
                          [Out]array<float>^% r, [Out]array<float>^% c,
                          [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax);
    static __int64 dgbequ(LapackLayout Layout, int m, int n, int kl, int ku,
                          array<double>^ ab, int ldab,
                          [Out]array<double>^% r, [Out]array<double>^% c,
                          [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax);

    static __int64 sgbequb(LapackLayout Layout, int m, int n, int kl, int ku,
                           array<float>^ ab, int ldab,
                           [Out]array<float>^% r, [Out]array<float>^% c,
                           [Out]float% rowCnd, [Out]float% colCnd, [Out]float% aMax);
    static __int64 dgbequb(LapackLayout Layout, int m, int n, int kl, int ku,
                           array<double>^ ab, int ldab,
                           [Out]array<double>^% r, [Out]array<double>^% c,
                           [Out]double% rowCnd, [Out]double% colCnd, [Out]double% aMax);

    static __int64 sgbcon(LapackLayout Layout, LapackNorm Norm, int n, int kl, int ku,
                          array<float>^ ab, int ldab, array<__int64>^ ipiv,
                          float aNorm, [Out]float% rCond);
    static __int64 dgbcon(LapackLayout Layout, LapackNorm Norm, int n, int kl, int ku,
                          array<double>^ ab, int ldab, array<__int64>^ ipiv,
                          double aNorm, [Out]double% rCond);

    static __int64 sgbrfs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int kl, int ku, int nrhs, array<float>^ ab, int ldab,
                          array<float>^ afb, int ldafb, array<__int64>^ ipiv,
                          array<float>^ b, int ldb, 
                          array<float>^ x, int ldx, 
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dgbrfs(LapackLayout Layout, LapackTranspose Trans,
                          int n, int kl, int ku, int nrhs, array<double>^ ab, int ldab,
                          array<double>^ afb, int ldafb, array<__int64>^ ipiv,
                          array<double>^ b, int ldb, 
                          array<double>^ x, int ldx, 
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 sgbrfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                           int n, int kl, int ku, int nrhs, array<float>^ ab, int ldab,
                           array<float>^ afb, int ldafb, array<__int64>^ ipiv,
                           array<float>^ r, array<float>^ c,
                           array<float>^ b, int ldb, 
                           array<float>^ x, int ldx, 
                           [Out]float% rCond, [Out]array<float>^% bErr,
                           int nErrBnds, [Out]array<float>^% errBndsNorm,
                           [Out]array<float>^% errBndsComp, int nParams, array<float>^ params);
    static __int64 dgbrfsx(LapackLayout Layout, LapackTranspose Trans, LapackEquil Equed,
                           int n, int kl, int ku, int nrhs, array<double>^ ab, int ldab,
                           array<double>^ afb, int ldafb, array<__int64>^ ipiv,
                           array<double>^ r, array<double>^ c,
                           array<double>^ b, int ldb, 
                           array<double>^ x, int ldx, 
                           [Out]double% rCond, [Out]array<double>^% bErr,
                           int nErrBnds, [Out]array<double>^% errBndsNorm,
                           [Out]array<double>^% errBndsComp, int nParams, array<double>^ params);
  #pragma endregion
  #pragma region general tridiagonal
    static __int64 sgttrf(int n, array<float>^ dl, array<float>^ d, array<float>^ du,
                          [Out]array<float>^% du2, [Out]array<__int64>^% ipiv);
    static __int64 dgttrf(int n, array<double>^ dl, array<double>^ d, array<double>^ du,
                          [Out]array<double>^% du2, [Out]array<__int64>^% ipiv);

    static __int64 sgttrs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                          array<float>^ dl, array<float>^ d, array<float>^ du,
                          array<float>^ du2, array<__int64>^ ipiv,
                          array<float>^ b, int ldb);
    static __int64 dgttrs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                          array<double>^ dl, array<double>^ d, array<double>^ du,
                          array<double>^ du2, array<__int64>^ ipiv,
                          array<double>^ b, int ldb);

    static __int64 sgtcon(LapackNorm Norm, int n,
                          array<float>^ dl, array<float>^ d, array<float>^ du,
                          array<float>^ du2, array<__int64>^ ipiv,
                          float aNorm, [Out]float% rCond);
    static __int64 dgtcon(LapackNorm Norm, int n,
                          array<double>^ dl, array<double>^ d, array<double>^ du,
                          array<double>^ du2, array<__int64>^ ipiv,
                          double aNorm, [Out]double% rCond);

    static __int64 sgtrfs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                          array<float>^ dl, array<float>^ d, array<float>^ du,
                          array<float>^ dlf, array<float>^ df, array<float>^ duf,
                          array<float>^ du2, array<__int64>^ ipiv,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dgtrfs(LapackLayout Layout, LapackTranspose Trans, int n, int nrhs,
                          array<double>^ dl, array<double>^ d, array<double>^ du,
                          array<double>^ dlf, array<double>^ df, array<double>^ duf,
                          array<double>^ du2, array<__int64>^ ipiv,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);
  #pragma endregion
  #pragma region diagonally dominant tridiagonal
    static __int64 sdttrfb(int n, array<float>^ dl, array<float>^ d, array<float>^ du);
    static __int64 ddttrfb(int n, array<double>^ dl, array<double>^ d, array<double>^ du);

    static __int64 sdttrsb(LapackTranspose Trans, int n, int nrhs,
                           array<float>^ dl, array<float>^ d, array<float>^ du,
                           array<float>^ b, int ldb);
    static __int64 ddttrsb(LapackTranspose Trans, int n, int nrhs,
                           array<double>^ dl, array<double>^ d, array<double>^ du,
                           array<double>^ b, int ldb);
  #pragma endregion
  #pragma region symmetric positive-definite
    static __int64 spotrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda);
    static __int64 dpotrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda);

    static __int64 spotrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ b, int ldb);
    static __int64 dpotrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ b, int ldb);

    static __int64 spoequ(LapackLayout Layout, int n,
                          array<float>^ a, int lda,
                          [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax);
    static __int64 dpoequ(LapackLayout Layout, int n,
                          array<double>^ a, int lda,
                          [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax);

    static __int64 spoequb(LapackLayout Layout, int n,
                           array<float>^ a, int lda,
                           [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax);
    static __int64 dpoequb(LapackLayout Layout, int n,
                           array<double>^ a, int lda,
                           [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax);

    static __int64 spocon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda,
                          float aNorm, [Out]float% rCond);
    static __int64 dpocon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda,
                          double aNorm, [Out]double% rCond);

    static __int64 sporfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dporfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 sporfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                           int n, int nrhs, array<float>^ a, int lda,
                           array<float>^ af, int ldaf, array<float>^ s,
                           array<float>^ b, int ldb,
                           array<float>^ x, int ldx,
                           [Out]float% rCond, [Out]array<float>^% bErr,
                           int nErrBnds, [Out]array<float>^% errBndsNorm,
                           [Out]array<float>^% errBndsComp,
                           int nParams, array<float>^ params);
    static __int64 dporfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                           int n, int nrhs, array<double>^ a, int lda,
                           array<double>^ af, int ldaf, array<double>^ s,
                           array<double>^ b, int ldb,
                           array<double>^ x, int ldx,
                           [Out]double% rCond, [Out]array<double>^% bErr,
                           int nErrBnds, [Out]array<double>^% errBndsNorm,
                           [Out]array<double>^% errBndsComp,
                           int nParams, array<double>^ params);
  #pragma endregion
  };
}