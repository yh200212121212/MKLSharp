#pragma once

using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Blas2 {
  public:
    static void sgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, long kl, long ku,
                      float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX, float beta, array<float>^ y, long incY);
    static void dgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, long kl, long ku,
                      double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX, double beta, array<double>^ y, long incY);

    static void sgemv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX, float beta, array<float>^ y, long incY);
    static void dgemv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX, double beta, array<double>^ y, long incY);

    static void sger(CBlasLayout Layout, long m, long n,
                     float alpha, array<float>^ x, long incX,
                     array<float>^ y, long incY,
                     array<float>^ a, long lda);
    static void dger(CBlasLayout Layout, long m, long n,
                     double alpha, array<double>^ x, long incX,
                     array<double>^ y, long incY,
                     array<double>^ a, long lda);

    static void ssbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, long k, float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX,
                      float beta, array<float>^ y, long incY);
    static void dsbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, long k, double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX,
                      double beta, array<double>^ y, long incY);

    static void sspmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, float alpha, array<float>^ ap,
                      array<float>^ x, long incX,
                      float beta, array<float>^ y, long incY);
    static void dspmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, double alpha, array<double>^ ap,
                      array<double>^ x, long incX,
                      double beta, array<double>^ y, long incY);

    static void sspr(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, float alpha, array<float>^ x, long incX, array<float>^ ap);
    static void dspr(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, double alpha, array<double>^ x, long incX, array<double>^ ap);

    static void sspr2(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, float alpha, array<float>^ x, long incX,
                     array<float>^ y, long incY, array<float>^ ap);
    static void dspr2(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, double alpha, array<double>^ x, long incX,
                     array<double>^ y, long incY, array<double>^ ap);

    static void ssymv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX,
                      float beta, array<float>^ y, long incY);
    static void dsymv(CBlasLayout Layout, CBlasUpLo UpLo,
                      long n, double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX,
                      double beta, array<double>^ y, long incY);

    static void ssyr(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, float alpha, array<float>^ x, long incX,
                     array<float>^ a, long lda);
    static void dsyr(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, double alpha, array<double>^ x, long incX,
                     array<double>^ a, long lda);

    static void ssyr2(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, float alpha, array<float>^ x, long incX,
                     array<float>^ y, long incY,
                     array<float>^ a, long lda);
    static void dsyr2(CBlasLayout Layout, CBlasUpLo UpLo,
                     long n, double alpha, array<double>^ x, long incX,
                     array<double>^ y, long incY,
                     array<double>^ a, long lda);

    static void stbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, long k, array<float>^ a, long lda,
                      array<float>^ x, long incX);
    static void dtbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, long k, array<double>^ a, long lda,
                      array<double>^ x, long incX);

    static void stbsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, long k, array<float>^ a, long lda,
                      array<float>^ x, long incX);
    static void dtbsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, long k, array<double>^ a, long lda,
                      array<double>^ x, long incX);

    static void stpmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<float>^ ap,
                      array<float>^ x, long incX);
    static void dtpmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<double>^ ap,
                      array<double>^ x, long incX);

    static void stpsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<float>^ ap,
                      array<float>^ x, long incX);
    static void dtpsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<double>^ ap,
                      array<double>^ x, long incX);

    static void strmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<float>^ a, long lda,
                      array<float>^ x, long incX);
    static void dtrmv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<double>^ a, long lda,
                      array<double>^ x, long incX);

    static void strsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<float>^ a, long lda,
                      array<float>^ x, long incX);
    static void dtrsv(CBlasLayout Layout, CBlasUpLo UpLo,
                      CBlasTranspose Trans, CBlasDiag Diag,
                      long n, array<double>^ a, long lda,
                      array<double>^ x, long incX);
  };
}