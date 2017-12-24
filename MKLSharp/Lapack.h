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

    static __int64 spotri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda);
    static __int64 dpotri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda);
  #pragma endregion
  #pragma region symmetric positive-definite, packed storage
    static __int64 spptrf(LapackLayout Layout, LapackUpLo UpLo, int n, array<float>^ ap);
    static __int64 dpptrf(LapackLayout Layout, LapackUpLo UpLo, int n, array<double>^ ap);

    static __int64 spptrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ ap, array<float>^ b, int ldb);
    static __int64 dpptrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ ap, array<double>^ b, int ldb);

    static __int64 sppequ(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ ap,
                          [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax);
    static __int64 dppequ(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ ap,
                          [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax);

    static __int64 sppcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ ap, float aNorm, [Out]float% rCond);
    static __int64 dppcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ ap, double aNorm, [Out]double% rCond);

    static __int64 spprfs(LapackLayout Layout, LapackUpLo UpLo, 
                          int n, int nrhs, array<float>^ ap, array<float>^ afp,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dpprfs(LapackLayout Layout, LapackUpLo UpLo, 
                          int n, int nrhs, array<double>^ ap, array<double>^ afp,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 spptri(LapackLayout Layout, LapackUpLo UpLo, int n, array<float>^ ap);
    static __int64 dpptri(LapackLayout Layout, LapackUpLo UpLo, int n, array<double>^ ap);
  #pragma endregion
  #pragma region symmetric positive-definite, RFP storage
    static __int64 spftrf(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, array<float>^ a);
    static __int64 dpftrf(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, array<double>^ a);

    static __int64 spftrs(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ a, array<float>^ b, int ldb);
    static __int64 dpftrs(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ a, array<double>^ b, int ldb);

    static __int64 spftri(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, array<float>^ a);
    static __int64 dpftri(LapackLayout Layout, LapackTranspose Trans, LapackUpLo UpLo,
                          int n, array<double>^ a);
  #pragma endregion
  #pragma region symmetric positive-definite, band
    static __int64 spbtrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<float>^ ab, int ldab);
    static __int64 dpbtrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<double>^ ab, int ldab);

    static __int64 spbtrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, int nrhs, array<float>^ ab, int ldab,
                          array<float>^ b, int ldb);
    static __int64 dpbtrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, int nrhs, array<double>^ ab, int ldab,
                          array<double>^ b, int ldb);

    static __int64 spbequ(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<float>^ ab, int ldab,
                          [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax);
    static __int64 dpbequ(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<double>^ ab, int ldab,
                          [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax);

    static __int64 spbcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<float>^ ab, int ldab,
                          float aNorm, [Out]float% rCond);
    static __int64 dpbcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, array<double>^ ab, int ldab,
                          double aNorm, [Out]double% rCond);

    static __int64 spbrfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, int nrhs, array<float>^ ab, int ldab,
                          array<float>^ afb, int ldafb, 
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dpbrfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int kd, int nrhs, array<double>^ ab, int ldab,
                          array<double>^ afb, int ldafb, 
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);
  #pragma endregion
  #pragma region symmetric positive-definite, tridiagonal
    static __int64 spttrf(int n, array<float>^ d, array<float>^ e);
    static __int64 dpttrf(int n, array<double>^ d, array<double>^ e);

    static __int64 spttrs(LapackLayout Layout, int n, int nrhs,
                          array<float>^ d, array<float>^ e,
                          array<float>^ b, int ldb);
    static __int64 dpttrs(LapackLayout Layout, int n, int nrhs,
                          array<double>^ d, array<double>^ e,
                          array<double>^ b, int ldb);

    static __int64 sptcon(int n, array<float>^ d, array<float>^ e,
                          float aNorm, [Out]float% rCond);
    static __int64 dptcon(int n, array<double>^ d, array<double>^ e,
                          double aNorm, [Out]double% rCond);

    static __int64 sptrfs(LapackLayout Layout, int n, int nrhs,
                          array<float>^ d, array<float>^ e,
                          array<float>^ df, array<float>^ ef,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dptrfs(LapackLayout Layout, int n, int nrhs,
                          array<double>^ d, array<double>^ e,
                          array<double>^ df, array<double>^ ef,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);
  #pragma endregion
  #pragma region symmetric indefinite
    static __int64 ssytrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda, [Out]array<__int64>^ ipiv);
    static __int64 dsytrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda, [Out]array<__int64>^ ipiv);

    static __int64 ssytrf_aa(LapackLayout Layout, LapackUpLo UpLo,
                             int n, array<float>^ a, int lda, [Out]array<__int64>^ ipiv);
    static __int64 dsytrf_aa(LapackLayout Layout, LapackUpLo UpLo,
                             int n, array<double>^ a, int lda, [Out]array<__int64>^ ipiv);

    static __int64 ssytrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                          array<float>^ b, int ldb);
    static __int64 dsytrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                          array<double>^ b, int ldb);

    static __int64 ssytrs2(LapackLayout Layout, LapackUpLo UpLo,
                           int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                           array<float>^ b, int ldb);
    static __int64 dsytrs2(LapackLayout Layout, LapackUpLo UpLo,
                           int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                           array<double>^ b, int ldb);

    static __int64 ssytrs_aa(LapackLayout Layout, LapackUpLo UpLo,
                             int n, int nrhs, array<float>^ a, int lda, array<__int64>^ ipiv,
                             array<float>^ b, int ldb);
    static __int64 dsytrs_aa(LapackLayout Layout, LapackUpLo UpLo,
                             int n, int nrhs, array<double>^ a, int lda, array<__int64>^ ipiv,
                             array<double>^ b, int ldb);

    static __int64 ssyequb(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<float>^ a, int lda,
                           [Out]array<float>^% s, [Out]float% sCond, [Out]float% aMax);
    static __int64 dsyequb(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<double>^ a, int lda,
                           [Out]array<double>^% s, [Out]double% sCond, [Out]double% aMax);

    static __int64 ssycon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda, array<__int64>^ ipiv,
                          float aNorm, [Out]float% rCond);
    static __int64 dsycon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda, array<__int64>^ ipiv,
                          double aNorm, [Out]double% rCond);

    static __int64 ssyrfs(LapackLayout Layout, LapackUpLo UpLo, 
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ af, int ldaf, array<__int64>^ ipiv,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx, 
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dsyrfs(LapackLayout Layout, LapackUpLo UpLo, 
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ af, int ldaf, array<__int64>^ ipiv,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx, 
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 ssyrfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                           int n, int nrhs, array<float>^ a, int lda,
                           array<float>^ af, int ldaf, array<__int64>^ ipiv,
                           array<float>^ s, array<float>^ b, int ldb,
                           array<float>^ x, int ldx,
                           [Out]float% rCond, [Out]array<float>^% bErr, 
                           int nErrBnds, [Out]array<float>^% errBndsNorm, 
                           [Out]array<float>^% errBndsComp,
                           int nParams, array<float>^ params);
    static __int64 dsyrfsx(LapackLayout Layout, LapackUpLo UpLo, LapackEquil Equed,
                           int n, int nrhs, array<double>^ a, int lda,
                           array<double>^ af, int ldaf, array<__int64>^ ipiv,
                           array<double>^ s, array<double>^ b, int ldb,
                           array<double>^ x, int ldx,
                           [Out]double% rCond, [Out]array<double>^% bErr, 
                           int nErrBnds, [Out]array<double>^% errBndsNorm, 
                           [Out]array<double>^% errBndsComp,
                           int nParams, array<double>^ params);

    static __int64 ssytri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ a, int lda, array<__int64>^ ipiv);
    static __int64 dsytri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ a, int lda, array<__int64>^ ipiv);

    static __int64 ssytri2(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<float>^ a, int lda, array<__int64>^ ipiv);
    static __int64 dsytri2(LapackLayout Layout, LapackUpLo UpLo,
                           int n, array<double>^ a, int lda, array<__int64>^ ipiv);

    static __int64 ssytri2x(LapackLayout Layout, LapackUpLo UpLo,
                            int n, array<float>^ a, int lda, array<__int64>^ ipiv, int nb);
    static __int64 dsytri2x(LapackLayout Layout, LapackUpLo UpLo,
                            int n, array<double>^ a, int lda, array<__int64>^ ipiv, int nb);
  #pragma endregion
  #pragma region symmetric indefinite, packed storage
    static __int64 ssptrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ ap, [Out]array<__int64>^ ipiv);
    static __int64 dsptrf(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ ap, [Out]array<__int64>^ ipiv);

    static void sspffrt2(int n, int nColumn, array<float>^ ap);
    static void dspffrt2(int n, int nColumn, array<double>^ ap);

    static void sspffrtx(int n, int nColumn, array<float>^ ap);
    static void dspffrtx(int n, int nColumn, array<double>^ ap);

    static __int64 ssptrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ ap, array<__int64>^ ipiv,
                          array<float>^ b, int ldb);
    static __int64 dsptrs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ ap, array<__int64>^ ipiv,
                          array<double>^ b, int ldb);

    static __int64 sspcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ ap, array<__int64>^ ipiv,
                          float aNorm, [Out]float% rCond);
    static __int64 dspcon(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ ap, array<__int64>^ ipiv,
                          double aNorm, [Out]double% rCond);

    static __int64 ssprfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<float>^ ap,
                          array<float>^ afp, array<__int64>^ ipiv,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx,
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dsprfs(LapackLayout Layout, LapackUpLo UpLo,
                          int n, int nrhs, array<double>^ ap,
                          array<double>^ afp, array<__int64>^ ipiv,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx,
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 ssptri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<float>^ ap, array<__int64>^ ipiv);
    static __int64 dsptri(LapackLayout Layout, LapackUpLo UpLo,
                          int n, array<double>^ ap, array<__int64>^ ipiv);
  #pragma endregion
  #pragma region triangular
    static __int64 strtrs(LapackLayout Layout, LapackUpLo UpLo,
                          LapackTranspose Trans, LapackDiag Diag,
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ b, int ldb);
    static __int64 dtrtrs(LapackLayout Layout, LapackUpLo UpLo,
                          LapackTranspose Trans, LapackDiag Diag,
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ b, int ldb);

    static __int64 strcon(LapackLayout Layout, LapackNorm Norm,
                          LapackUpLo UpLo, LapackDiag Diag, 
                          int n, array<float>^ a, int lda, [Out]float% rCond);
    static __int64 dtrcon(LapackLayout Layout, LapackNorm Norm,
                          LapackUpLo UpLo, LapackDiag Diag, 
                          int n, array<double>^ a, int lda, [Out]double% rCond);

    static __int64 strrfs(LapackLayout Layout, LapackUpLo UpLo,
                          LapackTranspose Trans, LapackDiag Diag, 
                          int n, int nrhs, array<float>^ a, int lda,
                          array<float>^ b, int ldb,
                          array<float>^ x, int ldx, 
                          [Out]array<float>^% fErr, [Out]array<float>^% bErr);
    static __int64 dtrrfs(LapackLayout Layout, LapackUpLo UpLo,
                          LapackTranspose Trans, LapackDiag Diag, 
                          int n, int nrhs, array<double>^ a, int lda,
                          array<double>^ b, int ldb,
                          array<double>^ x, int ldx, 
                          [Out]array<double>^% fErr, [Out]array<double>^% bErr);

    static __int64 strtri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                          int n, array<float>^ a, int lda);
    static __int64 dtrtri(LapackLayout Layout, LapackUpLo UpLo, LapackDiag Diag,
                          int n, array<double>^ a, int lda);
  #pragma endregion
  };
}